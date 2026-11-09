#pragma once
#include <filesystem>
#include <toml++/toml.hpp>

#include "flightsim/analysis/trim.hpp"

namespace flightsim::sim {
struct AircraftScenario {
  F16Model model;
  TrimSpec trim;
};
F16Model load_aircraft(const toml::table&, const std::filesystem::path& base_dir);
TrimSpec load_trim_spec(const toml::table&, double xcg);
AircraftScenario load_aircraft_scenario(const std::filesystem::path&);
}  // namespace flightsim::sim
