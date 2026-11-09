#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "flightsim/control/metrics.hpp"
#include "flightsim/sim/aircraft_scenario.hpp"
#include "flightsim/sim/runner.hpp"

namespace flightsim::sim {
// Closed-loop telemetry at each recorded sample (SI): the setpoints in effect and the cascade intermediates.
struct ApTrace {
  double h_cmd, V_cmd, theta_cmd, q_cmd;
};
// Peaks over EVERY integration step (not only the recorded samples). Elevator values come from the actuator
// state (all zero movement for ideal surfaces); throttle is the command.
struct ActuatorUsage {
  double de_position_peak_rad = 0, de_rate_peak_radps = 0, throttle_min = 1e300, throttle_max = -1e300;
};
struct AircraftRunOutput {
  std::vector<double> t;
  std::vector<Vec<17>> x;
  std::vector<F16Input>
      u;  // input applied over the step that STARTS at t (last sample repeats the final one)
  StopReason reason = StopReason::FinalTime;
  double t_stop = 0;
  bool out_of_domain = false;  // any table lookup left its domain (extrapolated); see first_violation
  std::string first_violation, message;
  F16Input u_trim;          // base controls (trim or explicit)
  std::vector<ApTrace> ap;  // aligned with t; empty unless the scenario has an [autopilot]
  ActuatorUsage usage;
};

// Control-performance report of a closed-loop run (empty without an [autopilot]); all SI. Computed from the
// recorded samples by the pure functions in control/metrics.hpp.
struct StepReport {
  std::string kind;  // "altitude" | "airspeed"
  double t_cmd_s, from, to;
  StepMetrics metrics;
  double peak_cross_deviation;  // the other loop's peak |error| (airspeed for an altitude step, altitude for
                                // an airspeed step)
};
struct DisturbanceReport {
  std::string kind;     // "wind_step" | "gust"
  double t_s, t_end_s;  // window: from this disturbance to the next one (or the end of the run)
  DisturbanceMetrics altitude, airspeed;  // deviation from the commanded altitude / airspeed
};
struct ControlReport {
  std::vector<StepReport> steps;
  std::vector<DisturbanceReport> disturbances;
};
ControlReport control_report(const AircraftScenario&, const AircraftRunOutput&);

// Open-loop fixed-step run of an aircraft scenario with a [sim] table. Throws if from_trim and trim fails.
// Stops at t_final, at the first boundary with h <= 0 (stop_at_ground), or on a non-finite / out-of-envelope
// state (NonFinite, with the cause in `message`).
AircraftRunOutput run_aircraft(const AircraftScenario&);

// Columns: t_s, x17 in SI, then V_mps..h_m (x13 chart), mach, qbar_Pa, thrust_N, throttle, de/da/dr_rad,
// then S&L-unit copies (V_fps, alpha_deg, ..., h_ft, pN_ft, pE_ft, de/da/dr_deg) named like the JSBSim CSVs.
void write_aircraft_csv(const std::filesystem::path&, const AircraftScenario&, const AircraftRunOutput&);
void write_aircraft_metadata(const std::filesystem::path&, const AircraftScenario&, const AircraftRunOutput&);
}  // namespace flightsim::sim
