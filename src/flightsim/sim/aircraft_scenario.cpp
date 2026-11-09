#include "flightsim/sim/aircraft_scenario.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include "flightsim/math/units.hpp"

namespace flightsim::sim {
namespace {
void keys(const toml::table& t, std::initializer_list<std::string_view> allowed) {
  for (const auto& [k, v] : t) {
    if (std::find(allowed.begin(), allowed.end(), k.str()) == allowed.end()) {
      throw std::runtime_error("unknown aircraft scenario key: " + std::string(k.str()));
    }
  }
}
const toml::table& section(const toml::table& t, const char* k) {
  const auto* p = t[k].as_table();
  if (!p) {
    throw std::runtime_error(std::string(k) + ": expected table");
  }
  return *p;
}
template <class T>
T required(const toml::table& t, const char* k) {
  const auto v = t[k].value<T>();
  if (!v) {
    throw std::runtime_error(std::string(k) + ": missing or invalid value");
  }
  if constexpr (std::is_same_v<T, double>) {
    if (!std::isfinite(*v)) {
      throw std::runtime_error(std::string(k) + ": must be finite");
    }
  }
  return *v;
}
Vec3 vec3(const toml::table& t, const char* k) {
  const auto* a = t[k].as_array();
  if (!a || a->size() != 3) {
    throw std::runtime_error(std::string(k) + ": expected an array of 3 numbers");
  }
  Vec3 v;
  for (std::size_t i = 0; i < 3; ++i) {
    const auto x = a->get(i)->value<double>();
    if (!x || !std::isfinite(*x)) {
      throw std::runtime_error(std::string(k) + ": expected finite numbers");
    }
    v(static_cast<Eigen::Index>(i)) = *x;
  }
  return v;
}
// Event times must sit on step boundaries (PRD 9.4): a held-constant-per-step input cannot represent others.
void require_aligned(double x, double dt, const char* what) {
  if (std::abs(x / dt - std::round(x / dt)) > 1e-9) {
    throw std::runtime_error(std::string(what) + ": must be a multiple of sim.dt_s");
  }
}
double dimension(const toml::table& t, const char* si, const char* imperial) {
  if (t.contains(si) == t.contains(imperial)) {
    throw std::runtime_error(std::string(si) + " or " + imperial + ": specify exactly one");
  }
  return t.contains(si) ? required<double>(t, si) : required<double>(t, imperial) * units::ft;
}
// A relative path is looked up next to the scenario, then in each parent directory (repo-root-relative
// paths).
std::filesystem::path resolve_path(const std::filesystem::path& base, std::filesystem::path path) {
  if (path.is_relative()) {
    auto candidate = base / path;
    if (!std::filesystem::exists(candidate)) {
      for (auto p = base; !p.empty(); p = p.parent_path()) {
        if (std::filesystem::exists(p / path)) {
          candidate = p / path;
          break;
        }
        if (p == p.root_path()) {
          break;
        }
      }
    }
    path = candidate;
  }
  return path;
}
}  // namespace
F16Model load_aircraft(const toml::table& t, const std::filesystem::path& base) {
  keys(t, {"aircraft", "trim", "initial", "inputs", "wind", "sim", "autopilot", "commands"});
  const auto& a = section(t, "aircraft");
  keys(a, {"data", "xcg", "atmosphere", "gravity", "rotor", "actuators"});
  const std::filesystem::path path = resolve_path(base, required<std::string>(a, "data"));
  F16Model m{load_f16_data(path)};
  m.xcg = required<double>(a, "xcg");
  const auto atm = required<std::string>(a, "atmosphere");
  if (atm == "sl_adc") {
    m.atm = AtmosphereKind::SlAdc;
  } else if (atm == "us76") {
    m.atm = AtmosphereKind::Us76;
  } else {
    throw std::runtime_error("atmosphere: expected sl_adc or us76");
  }
  if (const auto* g = a["gravity"].as_table()) {
    keys(*g, {"g0", "R"});
    m.grav.g0 = required<double>(*g, "g0");
    const auto radius = (*g)["R"].value<double>();
    if (!radius || std::isnan(*radius)) {
      throw std::runtime_error("gravity.R: expected positive radius or inf");
    }
    m.grav.R = *radius;
    if (m.grav.g0 <= 0 || m.grav.R <= 0) {
      throw std::runtime_error("gravity: g0/R must be positive");
    }
  } else {
    const auto kind = required<std::string>(a, "gravity");
    if (kind == "sl") {
      m.grav.g0 = m.data.g_sl;
    } else if (kind == "standard") {
      m.grav.g0 = units::g0;
    } else {
      throw std::runtime_error("gravity: expected sl, standard or {g0,R}");
    }
  }
  m.include_rotor = required<bool>(a, "rotor");
  m.actuators = required<bool>(a, "actuators");
  if (t.contains("wind")) {
    const auto& w = section(t, "wind");
    keys(w, {"constant_mps", "gust", "steps"});
    if (w.contains("constant_mps")) {
      m.wind.constant_N = vec3(w, "constant_mps");
    }
    if (w.contains("gust")) {
      const auto& g = section(w, "gust");
      keys(g, {"t_start_s", "duration_s", "amplitude_mps", "direction_N"});
      m.wind.gust = Gust{required<double>(g, "t_start_s"), required<double>(g, "duration_s"),
                         required<double>(g, "amplitude_mps"), vec3(g, "direction_N")};
    }
  }
  return m;
}
TrimSpec load_trim_spec(const toml::table& t, double cg) {
  keys(t, {"aircraft", "trim", "initial", "inputs", "wind", "sim", "autopilot", "commands"});
  const auto& r = section(t, "trim");
  keys(r, {"V_mps", "V_fps", "h_m", "h_ft", "gamma_deg"});
  TrimSpec s{dimension(r, "V_mps", "V_fps"), dimension(r, "h_m", "h_ft"),
             required<double>(r, "gamma_deg") * units::deg, cg};
  if (s.V <= 0) {
    throw std::runtime_error("trim.V: must be positive");
  }
  return s;
}
AircraftScenario load_aircraft_scenario(const std::filesystem::path& path) {
  const auto t = toml::parse_file(path.string());
  auto m = load_aircraft(t, std::filesystem::absolute(path).parent_path());
  auto s = load_trim_spec(t, m.xcg);
  AircraftScenario sc{std::move(m), s, {}, true, Vec<13>::Zero(), {}, {}, path.string(), {}, {}, {}};
  if (!t.contains("sim")) {
    if (t.contains("autopilot") || t.contains("commands")) {
      throw std::runtime_error("autopilot/commands: require a [sim] table");
    }
    return sc;
  }
  const auto& sm = section(t, "sim");
  keys(sm, {"integrator", "dt_s", "t_final_s", "output_every", "stop_at_ground"});
  AircraftRunConfig rc;
  rc.integ = integrator_from_string(required<std::string>(sm, "integrator"));
  rc.dt = required<double>(sm, "dt_s");
  rc.t_final = required<double>(sm, "t_final_s");
  rc.output_every =
      sm.contains("output_every") ? static_cast<int>(required<std::int64_t>(sm, "output_every")) : 1;
  rc.stop_at_ground = sm.contains("stop_at_ground") ? required<bool>(sm, "stop_at_ground") : true;
  if (rc.dt <= 0 || rc.t_final <= 0 || rc.output_every < 1) {
    throw std::runtime_error("sim: dt_s, t_final_s, output_every must be positive");
  }
  require_aligned(rc.t_final, rc.dt, "sim.t_final_s");
  if (const auto n = std::llround(rc.t_final / rc.dt); n % rc.output_every != 0) {
    throw std::runtime_error("sim.output_every: must divide t_final_s / dt_s");
  }
  sc.run = rc;
  if (const auto& g = sc.model.wind.gust) {
    require_aligned(g->t_start_s, rc.dt, "wind.gust.t_start_s");
    require_aligned(g->duration_s, rc.dt, "wind.gust.duration_s");
  }
  if (t.contains("initial")) {
    const auto& in = section(t, "initial");
    keys(in, {"from_trim", "x13", "throttle", "de_deg", "da_deg", "dr_deg"});
    sc.from_trim = required<bool>(in, "from_trim");
    if (!sc.from_trim) {
      const auto* a = in["x13"].as_array();
      if (!a || a->size() != 13) {
        throw std::runtime_error("initial.x13: expected 13 numbers (S&L order, SI)");
      }
      for (std::size_t i = 0; i < 13; ++i) {
        const auto x = a->get(i)->value<double>();
        if (!x || !std::isfinite(*x)) {
          throw std::runtime_error("initial.x13: expected finite numbers");
        }
        sc.x13_0(static_cast<Eigen::Index>(i)) = *x;
      }
      sc.u_base = {required<double>(in, "throttle"), required<double>(in, "de_deg") * units::deg,
                   required<double>(in, "da_deg") * units::deg, required<double>(in, "dr_deg") * units::deg};
    }
  }
  if (const auto* arr = t["inputs"].as_array()) {
    for (const auto& node : *arr) {
      const auto* e = node.as_table();
      if (!e) {
        throw std::runtime_error("inputs: expected an array of tables");
      }
      keys(*e, {"channel", "t_start_s", "t_end_s", "value", "mode"});
      if (e->contains("mode") && required<std::string>(*e, "mode") != "add") {
        throw std::runtime_error("inputs.mode: only 'add' is supported");
      }
      const auto ch = required<std::string>(*e, "channel");
      InputEvent ev{InputEvent::Throttle, required<double>(*e, "t_start_s"), required<double>(*e, "t_end_s"),
                    required<double>(*e, "value")};
      if (ch == "throttle") {
        ev.channel = InputEvent::Throttle;
      } else if (ch == "de" || ch == "da" || ch == "dr") {
        ev.channel = ch == "de" ? InputEvent::De : ch == "da" ? InputEvent::Da : InputEvent::Dr;
        ev.value *= units::deg;
      } else {
        throw std::runtime_error("inputs.channel: expected throttle, de, da or dr");
      }
      if (ev.t_end_s <= ev.t_start_s || ev.t_start_s < 0 || ev.t_end_s > rc.t_final) {
        throw std::runtime_error("inputs: need 0 <= t_start_s < t_end_s <= sim.t_final_s");
      }
      require_aligned(ev.t_start_s, rc.dt, "inputs.t_start_s");
      require_aligned(ev.t_end_s, rc.dt, "inputs.t_end_s");
      sc.inputs.push_back(ev);
    }
  }
  if (const auto* w = t["wind"].as_table(); w && w->contains("steps")) {
    const auto* arr = (*w)["steps"].as_array();
    if (!arr) {
      throw std::runtime_error("wind.steps: expected an array of tables");
    }
    for (const auto& node : *arr) {
      const auto* e = node.as_table();
      if (!e) {
        throw std::runtime_error("wind.steps: expected an array of tables");
      }
      keys(*e, {"t_s", "delta_mps"});
      const WindStep step{required<double>(*e, "t_s"), vec3(*e, "delta_mps")};
      if (step.t_s < 0 || step.t_s > rc.t_final) {
        throw std::runtime_error("wind.steps.t_s: must lie in [0, sim.t_final_s]");
      }
      require_aligned(step.t_s, rc.dt, "wind.steps.t_s");
      sc.wind_steps.push_back(step);
    }
  }
  if (t.contains("commands") && !t.contains("autopilot")) {
    throw std::runtime_error("commands: require an [autopilot] table");
  }
  if (t.contains("autopilot")) {
    if (!sc.inputs.empty()) {
      throw std::runtime_error("inputs: open-loop events cannot be combined with an [autopilot]");
    }
    const auto& ap = section(t, "autopilot");
    keys(ap, {"rate_hz", "gains", "antiwindup"});
    AutopilotConfig cfg;
    cfg.rate_hz = required<double>(ap, "rate_hz");
    cfg.antiwindup = required<bool>(ap, "antiwindup");
    const double ratio = 1 / (cfg.rate_hz * rc.dt);
    if (cfg.rate_hz <= 0 || ratio < 1 - 1e-9 || std::abs(ratio - std::round(ratio)) > 1e-9) {
      throw std::runtime_error("autopilot.rate_hz: sim.dt_s must divide 1/rate_hz (got dt_s=" +
                               std::to_string(rc.dt) + ", rate_hz=" + std::to_string(cfg.rate_hz) + ")");
    }
    cfg.steps_per_sample = std::lround(ratio);
    cfg.gains_path =
        resolve_path(std::filesystem::absolute(path).parent_path(), required<std::string>(ap, "gains"))
            .string();
    cfg.gains = load_autopilot_gains(cfg.gains_path);
    if (std::abs(cfg.gains.rate_hz - cfg.rate_hz) > 1e-9 * cfg.rate_hz) {
      throw std::runtime_error("autopilot.rate_hz: the gains were designed for " +
                               std::to_string(cfg.gains.rate_hz) + " Hz; redesign them for another rate");
    }
    sc.autopilot = cfg;
    if (const auto* arr = t["commands"].as_array()) {
      for (const auto& node : *arr) {
        const auto* e = node.as_table();
        if (!e) {
          throw std::runtime_error("commands: expected an array of tables");
        }
        keys(*e, {"kind", "t_s", "value_ft", "value_m", "value_fps", "value_mps"});
        const auto kind = required<std::string>(*e, "kind");
        Command c{Command::Altitude, required<double>(*e, "t_s"), 0};
        if (kind == "altitude") {
          c.value = dimension(*e, "value_m", "value_ft");
        } else if (kind == "airspeed") {
          c.kind = Command::Airspeed;
          c.value = dimension(*e, "value_mps", "value_fps");
        } else {
          throw std::runtime_error("commands.kind: expected altitude or airspeed");
        }
        if (c.t_s < 0 || c.t_s >= rc.t_final || c.value <= 0) {
          throw std::runtime_error("commands: need 0 <= t_s < sim.t_final_s and a positive value");
        }
        require_aligned(c.t_s, rc.dt, "commands.t_s");
        sc.commands.push_back(c);
      }
      std::stable_sort(sc.commands.begin(), sc.commands.end(),
                       [](const Command& a, const Command& b) { return a.t_s < b.t_s; });
    }
  }
  return sc;
}
Vec3 constant_wind_at(const AircraftScenario& s, double t) {
  Vec3 w = s.model.wind.constant_N;
  for (const auto& step : s.wind_steps) {
    if (t >= step.t_s - 1e-9) {  // event times are multiples of dt_s; slack only absorbs k*dt roundoff
      w += step.delta_N;
    }
  }
  return w;
}
}  // namespace flightsim::sim
