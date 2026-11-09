#include "flightsim/sim/aircraft_runner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <optional>
#include <stdexcept>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/control/autopilot.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/version.hpp"

#ifndef FLIGHTSIM_VERSION
#error "FLIGHTSIM_VERSION must be defined by the build"
#endif

namespace flightsim::sim {
namespace {
F16Input input_at(const AircraftScenario& s, const F16Input& base, double t) {
  F16Input u = base;
  for (const auto& e : s.inputs) {
    if (e.t_start_s <= t && t < e.t_end_s) {
      (e.channel == InputEvent::Throttle ? u.throttle
       : e.channel == InputEvent::De     ? u.de
       : e.channel == InputEvent::Da     ? u.da
                                         : u.dr) += e.value;
    }
  }
  return u;
}
}  // namespace

AircraftRunOutput run_aircraft(const AircraftScenario& s) {
  if (!s.run) {
    throw std::runtime_error("aircraft scenario has no [sim] table");
  }
  const auto& rc = *s.run;
  AircraftRunOutput out;
  Vec<17> x;
  if (s.from_trim) {
    const auto tr = trim_level(s.model, s.trim);
    if (!tr.ok) {
      throw std::runtime_error("trim failed: " + tr.message);
    }
    x = tr.x;
    out.u_trim = tr.u;
  } else {
    out.u_trim = s.u_base;
    x = x17_from_x13(s.model, s.x13_0, s.u_base);
  }
  LookupMonitor monitor;
  F16Model model = s.model;
  model.monitor = &monitor;

  const long n = std::lround(rc.t_final / rc.dt);
  const auto time_at = [&](long k) { return k == n ? rc.t_final : static_cast<double>(k) * rc.dt; };
  StepStats stats;

  // Closed loop: the controller is the one stateful component. It runs OUTSIDE every derivative evaluation,
  // once per sample before the integration step that starts on a sample instant, and its output is held.
  const auto measure = [&](double t, const Vec<17>& xs) {
    const Vec<13> c = x13_from_x17(model, xs, t);
    return ApMeasurement{c(11), c(0), c(4), c(7)};
  };
  std::optional<Autopilot> ap;
  ApCommand cmd{};
  ApOutput ap_out{};
  if (s.autopilot) {
    model.wind.constant_N = constant_wind_at(s, 0.0);
    const ApMeasurement m0 = measure(0.0, x);
    ap.emplace(s.autopilot->gains, s.autopilot->antiwindup);
    ap->engage(m0, out.u_trim);  // bumpless: the first sample outputs exactly the trim command
    cmd = {m0.h, m0.V};
    ap_out = {out.u_trim, m0.theta, 0.0};
  }
  const auto trace = [&] { return ApTrace{cmd.h, cmd.V, ap_out.theta_cmd, ap_out.q_cmd}; };
  const auto record = [&](double t, const Vec<17>& xs, const F16Input& uu) {
    out.t.push_back(t);
    out.x.push_back(xs);
    out.u.push_back(uu);
    if (ap) {
      out.ap.push_back(trace());
    }
  };
  F16Input u = ap ? out.u_trim : input_at(s, out.u_trim, 0.0);
  record(0.0, x, u);
  bool last_recorded = true;  // was the sample at the start of this step recorded?
  out.t_stop = rc.t_final;
  for (long k = 0; k < n; ++k) {
    const double t = time_at(k);
    model.wind.constant_N = constant_wind_at(s, t);
    const Vec<17> x_before = x;
    try {
      if (ap) {
        if (k % s.autopilot->steps_per_sample == 0) {
          for (const auto& c : s.commands) {  // sorted by time: the latest command of each kind wins
            if (k >= std::lround(c.t_s / rc.dt)) {
              (c.kind == Command::Altitude ? cmd.h : cmd.V) = c.value;
            }
          }
          ap_out = ap->update(measure(t, x), cmd);
        }
        u = ap_out.u;
        if (last_recorded) {  // the sample at t carries the command applied over the step that starts there
          out.u.back() = u;
          out.ap.back() = trace();
        }
      } else {
        u = input_at(s, out.u_trim, t);
      }
      x = step(
          rc.integ, [&](double tt, const Vec<17>& xx) { return model.deriv(tt, xx, u); }, t, x, rc.dt, true,
          stats);
    } catch (const std::exception& e) {  // atmosphere altitude domain, non-finite inputs
      out.reason = StopReason::NonFinite;
      out.t_stop = time_at(k + 1);
      out.message = std::string("step failed at t=") + std::to_string(t) + " s: " + e.what();
      break;
    }
    const double tn = time_at(k + 1);
    if (!x.allFinite()) {
      out.reason = StopReason::NonFinite;
      out.t_stop = tn;
      out.message = "non-finite state at t=" + std::to_string(tn) + " s";
      break;
    }
    auto& use = out.usage;
    use.de_position_peak_rad = std::max(use.de_position_peak_rad, std::abs(x(ac::kDe)));
    use.de_rate_peak_radps =
        std::max(use.de_rate_peak_radps, std::abs(x(ac::kDe) - x_before(ac::kDe)) / rc.dt);
    use.throttle_min = std::min(use.throttle_min, u.throttle);
    use.throttle_max = std::max(use.throttle_max, u.throttle);
    const bool ground = rc.stop_at_ground && -x(rb::kPos + 2) <= 0.0;
    last_recorded = ground || (k + 1) % rc.output_every == 0;
    if (last_recorded) {
      record(tn, x, ap ? u : input_at(s, out.u_trim, tn));
    }
    if (ground) {
      out.reason = StopReason::GroundContact;
      out.t_stop = tn;
      break;
    }
  }
  out.out_of_domain = monitor.out_of_domain;
  out.first_violation = monitor.first_violation;
  return out;
}

namespace {
constexpr const char* kCols[] = {
    "t_s",       "pN_m",         "pE_m",         "pD_m",         "vN_mps",  "vE_mps",    "vD_mps",
    "q0",        "q1",           "q2",           "q3",           "p_radps", "q_radps",   "r_radps",
    "power_pct", "de_state_rad", "da_state_rad", "dr_state_rad", "V_mps",   "alpha_rad", "beta_rad",
    "phi_rad",   "theta_rad",    "psi_rad",      "h_m",          "mach",    "qbar_Pa",   "thrust_N",
    "throttle",  "de_rad",       "da_rad",       "dr_rad",       "V_fps",   "alpha_deg", "beta_deg",
    "phi_deg",   "theta_deg",    "psi_deg",      "p_dps",        "q_dps",   "r_dps",     "h_ft",
    "pN_ft",     "pE_ft",        "de_deg",       "da_deg",       "dr_deg"};
}  // namespace

void write_aircraft_csv(const std::filesystem::path& path, const AircraftScenario& s,
                        const AircraftRunOutput& out) {
  std::ofstream f(path);
  if (!f) {
    throw std::runtime_error("cannot write " + path.string());
  }
  F16Model model = s.model;
  model.monitor = nullptr;
  bool first = true;
  for (const char* c : kCols) {
    f << (first ? "" : ",") << c;
    first = false;
  }
  if (!out.ap.empty()) {
    f << ",h_cmd_m,V_cmd_mps,theta_cmd_rad,q_cmd_radps";
  }
  f << '\n';
  char buf[40];
  for (std::size_t i = 0; i < out.t.size(); ++i) {
    const Vec<17>& x = out.x[i];
    const F16Input& u = out.u[i];
    model.wind.constant_N = constant_wind_at(s, out.t[i]);
    const auto ev = model.evaluate(out.t[i], x, u);
    const Vec<13> c = x13_from_x17(model, x, out.t[i]);
    const double d = units::deg, ft = units::ft;
    const double row[] = {out.t[i],   x(0),       x(1),     x(2),     x(3),      x(4),     x(5),
                          x(6),       x(7),       x(8),     x(9),     x(10),     x(11),    x(12),
                          x(13),      x(14),      x(15),    x(16),    c(0),      c(1),     c(2),
                          c(3),       c(4),       c(5),     c(11),    ev.mach,   ev.qbar,  ev.thrust,
                          u.throttle, u.de,       u.da,     u.dr,     c(0) / ft, c(1) / d, c(2) / d,
                          c(3) / d,   c(4) / d,   c(5) / d, c(6) / d, c(7) / d,  c(8) / d, c(11) / ft,
                          c(9) / ft,  c(10) / ft, u.de / d, u.da / d, u.dr / d};
    for (std::size_t j = 0; j < std::size(row); ++j) {
      std::snprintf(buf, sizeof buf, "%.17g", row[j]);
      f << (j ? "," : "") << buf;
    }
    if (!out.ap.empty()) {
      const ApTrace& a = out.ap[i];
      for (const double v : {a.h_cmd, a.V_cmd, a.theta_cmd, a.q_cmd}) {
        std::snprintf(buf, sizeof buf, "%.17g", v);
        f << ',' << buf;
      }
    }
    f << '\n';
  }
}

ControlReport control_report(const AircraftScenario& s, const AircraftRunOutput& out) {
  ControlReport r;
  if (!s.autopilot || !s.run || out.ap.size() != out.t.size()) {
    return r;
  }
  F16Model model = s.model;
  model.monitor = nullptr;
  std::vector<double> h, V;
  for (std::size_t i = 0; i < out.t.size(); ++i) {
    model.wind.constant_N = constant_wind_at(s, out.t[i]);
    const Vec<13> c = x13_from_x17(model, out.x[i], out.t[i]);
    h.push_back(c(11));
    V.push_back(c(0));
  }
  const double t_end = out.t.back();
  const auto first_at = [&](double time) {
    return static_cast<std::size_t>(std::lower_bound(out.t.begin(), out.t.end(), time - 1e-9) -
                                    out.t.begin());
  };
  // Each command's window runs to the next command (of either kind) or the end of the record.
  for (std::size_t j = 0; j < s.commands.size(); ++j) {
    const auto& c = s.commands[j];
    const double t_to = j + 1 < s.commands.size() ? s.commands[j + 1].t_s : t_end;
    const std::size_t i0 = first_at(c.t_s), i1 = first_at(t_to) + (j + 1 < s.commands.size() ? 0 : 1);
    const auto slice = [&](const std::vector<double>& v) {
      return std::vector<double>(v.begin() + static_cast<std::ptrdiff_t>(i0),
                                 v.begin() + static_cast<std::ptrdiff_t>(std::min(i1, v.size())));
    };
    const bool altitude = c.kind == Command::Altitude;
    // The setpoint before this command: the last earlier command of the same kind, else the trimmed value.
    double from = altitude ? out.ap.front().h_cmd : out.ap.front().V_cmd;
    for (std::size_t q = 0; q < j; ++q) {
      if (s.commands[q].kind == c.kind) {
        from = s.commands[q].value;
      }
    }
    const auto t = slice(out.t);
    StepReport sr{altitude ? "altitude" : "airspeed",
                  c.t_s,
                  from,
                  c.value,
                  step_metrics(t, slice(altitude ? h : V), c.t_s, from, c.value),
                  0.0};
    const double other_ref = altitude ? out.ap[i0].V_cmd : out.ap[i0].h_cmd;
    sr.peak_cross_deviation = peak_abs_deviation(t, slice(altitude ? V : h), other_ref, c.t_s, t_to);
    r.steps.push_back(sr);
  }
  // Disturbance windows run from one disturbance to the next (or the end of the record).
  std::vector<std::pair<double, std::string>> events;
  for (const auto& w : s.wind_steps) {
    events.emplace_back(w.t_s, "wind_step");
  }
  if (s.model.wind.gust) {
    events.emplace_back(s.model.wind.gust->t_start_s, "gust");
  }
  std::sort(events.begin(), events.end());
  for (std::size_t j = 0; j < events.size(); ++j) {
    const double t0 = events[j].first, t1 = j + 1 < events.size() ? events[j + 1].first : t_end;
    const ApTrace& ref = out.ap[first_at(t0)];
    r.disturbances.push_back({events[j].second, t0, t1, disturbance_metrics(out.t, h, ref.h_cmd, t0, t1),
                              disturbance_metrics(out.t, V, ref.V_cmd, t0, t1)});
  }
  return r;
}

namespace {
// JSON has no NaN: a metric that never settled/recovered is omitted and flagged with `<name>_reached =
// false`.
void put_time(toml::table& t, const char* name, double seconds) {
  t.insert(std::string(name) + "_reached", !std::isnan(seconds));
  if (!std::isnan(seconds)) {
    t.insert(name, seconds);
  }
}
toml::table disturbance_table(const DisturbanceMetrics& d) {
  toml::table t;
  t.insert("peak_deviation", d.peak_deviation);
  put_time(t, "recovery_2pct_s", d.recovery_2pct_s);
  return t;
}
toml::table control_metadata(const AircraftScenario& s, const AircraftRunOutput& out) {
  const auto rep = control_report(s, out);
  toml::table c;
  c.insert("rate_hz", s.autopilot->rate_hz);
  c.insert("antiwindup", s.autopilot->antiwindup);
  c.insert("gains", s.autopilot->gains_path);
  toml::array steps;
  for (const auto& r : rep.steps) {
    toml::table t;
    t.insert("kind", r.kind);
    t.insert("t_cmd_s", r.t_cmd_s);
    t.insert("from", r.from);
    t.insert("to", r.to);
    t.insert("step", r.metrics.step);
    t.insert("overshoot_pct", r.metrics.overshoot_pct);
    put_time(t, "settling_2pct_s", r.metrics.settling_2pct_s);
    t.insert("steady_state_error", r.metrics.steady_state_error);
    t.insert("peak_cross_deviation", r.peak_cross_deviation);
    steps.push_back(std::move(t));
  }
  c.insert("steps", std::move(steps));
  toml::array disturbances;
  for (const auto& d : rep.disturbances) {
    toml::table t;
    t.insert("kind", d.kind);
    t.insert("t_s", d.t_s);
    t.insert("t_end_s", d.t_end_s);
    t.insert("altitude", disturbance_table(d.altitude));
    t.insert("airspeed", disturbance_table(d.airspeed));
    disturbances.push_back(std::move(t));
  }
  c.insert("disturbances", std::move(disturbances));
  toml::table u;
  u.insert("elevator_position_peak_rad", out.usage.de_position_peak_rad);
  u.insert("elevator_rate_peak_radps", out.usage.de_rate_peak_radps);
  u.insert("throttle_min", out.usage.throttle_min);
  u.insert("throttle_max", out.usage.throttle_max);
  c.insert("actuator_usage", std::move(u));
  return c;
}
}  // namespace

void write_aircraft_metadata(const std::filesystem::path& path, const AircraftScenario& s,
                             const AircraftRunOutput& out) {
  const auto& rc = *s.run;
  toml::table m;
  m.insert("version", FLIGHTSIM_VERSION);
  m.insert("scenario", s.source_path);
  m.insert("integrator", to_string(rc.integ));
  m.insert("dt_s", rc.dt);
  m.insert("t_final_s", rc.t_final);
  m.insert("xcg", s.model.xcg);
  m.insert("stop_reason", to_string(out.reason));
  m.insert("t_stop_s", out.t_stop);
  m.insert("out_of_domain", out.out_of_domain);
  m.insert("first_violation", out.first_violation);
  m.insert("message", out.message);
  if (s.autopilot) {
    m.insert("control", control_metadata(s, out));
  }
  std::ofstream f(path);
  if (!f) {
    throw std::runtime_error("cannot write " + path.string());
  }
  f << toml::json_formatter{m} << '\n';
}
}  // namespace flightsim::sim
