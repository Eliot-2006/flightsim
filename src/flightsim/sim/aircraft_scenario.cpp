#include "flightsim/sim/aircraft_scenario.hpp"

#include <algorithm>
#include <cmath>
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
double dimension(const toml::table& t, const char* si, const char* imperial) {
  if (t.contains(si) == t.contains(imperial)) {
    throw std::runtime_error(std::string(si) + " or " + imperial + ": specify exactly one");
  }
  return t.contains(si) ? required<double>(t, si) : required<double>(t, imperial) * units::ft;
}
}  // namespace
F16Model load_aircraft(const toml::table& t, const std::filesystem::path& base) {
  keys(t, {"aircraft", "trim"});
  const auto& a = section(t, "aircraft");
  keys(a, {"data", "xcg", "atmosphere", "gravity", "rotor", "actuators"});
  std::filesystem::path path = required<std::string>(a, "data");
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
  return m;
}
TrimSpec load_trim_spec(const toml::table& t, double cg) {
  keys(t, {"aircraft", "trim"});
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
  return {std::move(m), s};
}
}  // namespace flightsim::sim
