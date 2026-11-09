#pragma once
#include <filesystem>
#include <optional>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

#include "flightsim/analysis/trim.hpp"
#include "flightsim/control/autopilot.hpp"
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
// [autopilot]: closed-loop engagement. The controller samples at rate_hz (dt_s divides 1/rate_hz, so a sample
// falls every `steps_per_sample` integration steps) and holds its outputs in between.
struct AutopilotConfig {
  double rate_hz = 0;
  long steps_per_sample = 0;
  bool antiwindup = true;
  AutopilotGains gains;
  std::string gains_path;
};
// [[commands]]: ABSOLUTE setpoint from t_s on, SI (altitude [m], air-relative speed [m/s]). Before the first
// command of a kind the setpoint is the trimmed value. Takes effect at the first controller sample at or
// after t_s.
struct Command {
  enum Kind { Altitude, Airspeed } kind;
  double t_s, value;
};
// [[wind.steps]]: `delta_N` [m/s] is ADDED to the constant wind from t_s on (a time-aligned constant-wind
// step; t_s is a multiple of dt_s, so no integration step straddles it).
struct WindStep {
  double t_s;
  Vec3 delta_N;
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
  std::optional<AutopilotConfig> autopilot;
  std::vector<Command> commands;  // sorted by t_s (stable)
  std::vector<WindStep> wind_steps;
};
// Constant wind in effect over the step that starts at t (base constant wind plus every step with t_s <= t).
Vec3 constant_wind_at(const AircraftScenario&, double t);
F16Model load_aircraft(const toml::table&, const std::filesystem::path& base_dir);
TrimSpec load_trim_spec(const toml::table&, double xcg);
AircraftScenario load_aircraft_scenario(const std::filesystem::path&);
}  // namespace flightsim::sim
