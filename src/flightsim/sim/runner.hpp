#pragma once

#include <string>
#include <vector>

#include "flightsim/sim/scenario.hpp"

namespace flightsim::sim {

enum class StopReason { FinalTime, GroundContact, NonFinite };

struct Sample {
  double t;
  Vec<13> x;
};

struct RunOutput {
  std::vector<Sample> samples;  // t = 0, every output_every steps, and the stop boundary; finite states only
  StopReason reason = StopReason::FinalTime;
  double t_stop = 0.0;  // time of the boundary where the run ended
  long nfev = 0;
  std::string message;  // set for NonFinite: time and offending state
};

// Fixed-step run. Stops at t_final, at the first boundary with h <= 0 (if stop_at_ground), or at the first
// boundary where the state is not finite (the offending state goes into `message`, not into `samples`).
RunOutput run(const RigidBodyScenario&);

const char* to_string(StopReason);

}  // namespace flightsim::sim
