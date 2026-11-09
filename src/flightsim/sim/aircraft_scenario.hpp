#pragma once
#include <filesystem>
#include <optional>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

#include "flightsim/analysis/trim.hpp"
#include "flightsim/core/integrators.hpp"

namespace flightsim::sim {
// One open-loop input event: `value` is ADDED to the base input of `channel` on t_start <= t < t_end.
// Stored in model units (throttle 0..1, surfaces rad); the TOML gives surfaces in degrees.
struct InputEvent {
  enum Channel { Throttle, De, Da, Dr } channel;
  double t_start_s, t_end_s, value;
};
struct AircraftRunConfig {
  Integrator integ = Integrator::RK4;
  double dt = 0, t_final = 0;
  int output_every = 1;
  bool stop_at_ground = true;
};
struct AircraftScenario {
  F16Model model;
  TrimSpec trim;
  // Present when the file has [sim]; trim-only scenarios (trim/linearize CLI) leave it empty.
  std::optional<AircraftRunConfig> run;
  bool from_trim = true;  // [initial] from_trim; else x13 (S&L order, SI) + base controls below
  Vec<13> x13_0 = Vec<13>::Zero();
  F16Input u_base;  // base controls when !from_trim
  std::vector<InputEvent> inputs;
  std::string source_path;
};
F16Model load_aircraft(const toml::table&, const std::filesystem::path& base_dir);
TrimSpec load_trim_spec(const toml::table&, double xcg);
AircraftScenario load_aircraft_scenario(const std::filesystem::path&);
}  // namespace flightsim::sim
