#pragma once
#include <filesystem>

#include "flightsim/aircraft/f16_model.hpp"
#include "flightsim/control/pid.hpp"

namespace flightsim {

// Gains exactly as written by python/control_design/design.py (scenarios/control/autopilot_gains.toml): SI
// units, angles in rad, PID outputs are ABSOLUTE commands (elevator rad, throttle 0..1, theta/q references).
struct AutopilotGains {
  double rate_hz = 0;
  PidGains pitch_rate{}, pitch{}, altitude{}, airspeed{};
  double airspeed_to_pitch = 0;       // k_V = d theta_trim / d V [rad/(m/s)]
  double reference_airspeed_mps = 0;  // V_trim the cross-feed is relative to
};
// Strict schema: unknown or missing keys throw std::runtime_error naming the key.
AutopilotGains load_autopilot_gains(const std::filesystem::path&);

struct ApMeasurement {
  double h, V, theta, q;  // altitude [m, up], AIR-relative speed [m/s], pitch [rad], body pitch rate [rad/s]
};
struct ApCommand {
  double h, V;  // [m], [m/s]
};
struct ApOutput {
  F16Input u;
  double theta_cmd, q_cmd;  // cascade intermediates (telemetry)
};

// The structure design.md documents (the one stateful component, called outside every derivative):
//   theta_cmd = clamp(PI_h(h_cmd - h) + k_V (V - V_ref), -15 deg, +15 deg)   altitude PI + airspeed
//   cross-feed q_cmd     = P(theta_cmd - theta) delta_e   = PI(q_cmd - q) (negative gain: S&L +de = nose
//   down) throttle  = PI(V_cmd - V)                                                 limited to [0,1]
// `update` is one controller sample of period 1/rate_hz; the caller holds the returned inputs until the next
// sample. The ±15 deg bound applies to the COMBINED theta_cmd, so the altitude PI's back-calculation tracks
// that final bound (not only its own).
class Autopilot {
 public:
  Autopilot(const AutopilotGains&, bool antiwindup);
  // Bumpless engage from trim: every PID is reset so that the next update outputs exactly the trim command
  // (altitude -> trim theta less the cross-feed, pitch -> q = 0, pitch rate -> trim elevator, airspeed ->
  // trim throttle). `m` is the trimmed measurement; aileron and rudder stay at u_trim's values.
  void engage(const ApMeasurement& m, const F16Input& u_trim);
  ApOutput update(const ApMeasurement& m, const ApCommand& cmd);  // throws std::logic_error before engage()

 private:
  AutopilotGains g_;
  Pid altitude_, pitch_, pitch_rate_, airspeed_;
  F16Input u_trim_{};
  bool engaged_ = false;
};
}  // namespace flightsim
