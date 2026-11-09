#include "flightsim/control/autopilot.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>
#include <toml++/toml.hpp>

namespace flightsim {
namespace {
void keys(const toml::table& t, std::initializer_list<std::string_view> allowed, const char* where) {
  for (const auto& [k, v] : t) {
    if (std::find(allowed.begin(), allowed.end(), k.str()) == allowed.end()) {
      throw std::runtime_error(std::string("unknown autopilot gains key: ") + where + "." +
                               std::string(k.str()));
    }
  }
}
double number(const toml::table& t, const char* k, const char* where) {
  const auto v = t[k].value<double>();
  if (!v || !std::isfinite(*v)) {
    throw std::runtime_error(std::string("autopilot gains: ") + where + "." + k + " missing or not finite");
  }
  return *v;
}
const toml::table& table(const toml::table& t, const char* k) {
  const auto* p = t[k].as_table();
  if (!p) {
    throw std::runtime_error(std::string("autopilot gains: missing table [") + k + "]");
  }
  return *p;
}
PidGains pid_gains(const toml::table& root, const char* name) {
  const auto& t = table(root, name);
  keys(t, {"kp", "ki", "kd", "tf", "tt", "umin", "umax"}, name);
  return {number(t, "kp", name), number(t, "ki", name),   number(t, "kd", name),  number(t, "tf", name),
          number(t, "tt", name), number(t, "umin", name), number(t, "umax", name)};
}
}  // namespace

AutopilotGains load_autopilot_gains(const std::filesystem::path& path) {
  toml::table root;
  try {
    root = toml::parse_file(path.string());
  } catch (const toml::parse_error& e) {
    throw std::runtime_error("autopilot gains " + path.string() + ": " + std::string(e.description()));
  }
  keys(root, {"rate_hz", "pitch_rate", "pitch", "altitude", "airspeed", "crossfeed"}, "");
  AutopilotGains g;
  g.rate_hz = number(root, "rate_hz", "");
  if (g.rate_hz <= 0) {
    throw std::runtime_error("autopilot gains: rate_hz must be positive");
  }
  g.pitch_rate = pid_gains(root, "pitch_rate");
  g.pitch = pid_gains(root, "pitch");
  g.altitude = pid_gains(root, "altitude");
  g.airspeed = pid_gains(root, "airspeed");
  const auto& cf = table(root, "crossfeed");
  keys(cf, {"airspeed_to_pitch", "reference_airspeed_mps"}, "crossfeed");
  g.airspeed_to_pitch = number(cf, "airspeed_to_pitch", "crossfeed");
  g.reference_airspeed_mps = number(cf, "reference_airspeed_mps", "crossfeed");
  return g;
}

Autopilot::Autopilot(const AutopilotGains& g, bool antiwindup)
    : g_(g), altitude_(g.altitude), pitch_(g.pitch), pitch_rate_(g.pitch_rate), airspeed_(g.airspeed) {
  for (Pid* p : {&altitude_, &pitch_, &pitch_rate_, &airspeed_}) {
    p->antiwindup = antiwindup;
  }
}

void Autopilot::engage(const ApMeasurement& m, const F16Input& u_trim) {
  u_trim_ = u_trim;
  // The altitude PI output is theta_cmd less the cross-feed, so that theta_cmd (not the PI) equals trim
  // theta.
  altitude_.reset(m.theta - g_.airspeed_to_pitch * (m.V - g_.reference_airspeed_mps));
  pitch_.reset(0);  // trim pitch rate
  pitch_rate_.reset(u_trim.de);
  airspeed_.reset(u_trim.throttle);
  engaged_ = true;
}

ApOutput Autopilot::update(const ApMeasurement& m, const ApCommand& cmd) {
  if (!engaged_) {
    throw std::logic_error("Autopilot::update before engage()");
  }
  const double T = 1 / g_.rate_hz;
  const double lo = g_.altitude.umin, hi = g_.altitude.umax;  // the +-15 deg bound on the COMBINED theta_cmd
  const double crossfeed = g_.airspeed_to_pitch * (m.V - g_.reference_airspeed_mps);
  // Limit the altitude PI so that PI + cross-feed stays inside [lo, hi]; its back-calculation then tracks the
  // final bound. (clamp is monotone, so lo_pi <= hi_pi always.)
  const double theta_pi =
      altitude_.update(cmd.h, m.h, T, std::clamp(lo - crossfeed, lo, hi), std::clamp(hi - crossfeed, lo, hi));
  ApOutput out;
  out.theta_cmd =
      std::clamp(theta_pi + crossfeed, lo, hi);  // exact design equation; guards roundoff/extremes
  out.q_cmd = pitch_.update(out.theta_cmd, m.theta, T);
  out.u = u_trim_;
  out.u.de = pitch_rate_.update(out.q_cmd, m.q, T);
  out.u.throttle = airspeed_.update(cmd.V, m.V, T);
  return out;
}
}  // namespace flightsim
