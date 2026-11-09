#pragma once
#include <filesystem>
#include <string>
#include <vector>

#include "flightsim/sim/aircraft_scenario.hpp"
#include "flightsim/sim/runner.hpp"

namespace flightsim::sim {
struct AircraftRunOutput {
  std::vector<double> t;
  std::vector<Vec<17>> x;
  std::vector<F16Input>
      u;  // input applied over the step that STARTS at t (last sample repeats the final one)
  StopReason reason = StopReason::FinalTime;
  double t_stop = 0;
  bool out_of_domain = false;  // any table lookup left its domain (extrapolated); see first_violation
  std::string first_violation, message;
  F16Input u_trim;  // base controls (trim or explicit)
};

// Open-loop fixed-step run of an aircraft scenario with a [sim] table. Throws if from_trim and trim fails.
// Stops at t_final, at the first boundary with h <= 0 (stop_at_ground), or on a non-finite / out-of-envelope
// state (NonFinite, with the cause in `message`).
AircraftRunOutput run_aircraft(const AircraftScenario&);

// Columns: t_s, x17 in SI, then V_mps..h_m (x13 chart), mach, qbar_Pa, thrust_N, throttle, de/da/dr_rad,
// then S&L-unit copies (V_fps, alpha_deg, ..., h_ft, pN_ft, pE_ft, de/da/dr_deg) named like the JSBSim CSVs.
void write_aircraft_csv(const std::filesystem::path&, const AircraftScenario&, const AircraftRunOutput&);
void write_aircraft_metadata(const std::filesystem::path&, const AircraftScenario&, const AircraftRunOutput&);
}  // namespace flightsim::sim
