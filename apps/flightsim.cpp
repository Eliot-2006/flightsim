#include <cstring>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "flightsim/sim/output.hpp"
#include "flightsim/sim/runner.hpp"
#include "flightsim/sim/scenario.hpp"
#include "flightsim/version.hpp"

namespace {

int usage() {
  std::cerr << "usage: flightsim --version\n"
               "       flightsim run <scenario.toml> --out <dir>\n";
  return 2;
}

// Writes <dir>/trajectory.csv and <dir>/metadata.json. Exit 0 on a clean stop, 1 if the state went
// non-finite.
int cmd_run(const std::string& scenario_path, const std::filesystem::path& out_dir) {
  namespace sim = flightsim::sim;
  const sim::RigidBodyScenario scenario = sim::load_rigid_body_scenario(scenario_path);
  const sim::RunOutput out = sim::run(scenario);
  std::filesystem::create_directories(out_dir);
  sim::write_csv(out_dir / "trajectory.csv", scenario, out);
  sim::write_metadata(out_dir / "metadata.json", scenario, out);
  std::cout << "wrote " << (out_dir / "trajectory.csv").string() << " and "
            << (out_dir / "metadata.json").string() << " (" << out.samples.size()
            << " samples, stop: " << sim::to_string(out.reason) << " at t=" << out.t_stop << " s)\n";
  if (out.reason == sim::StopReason::NonFinite) {
    std::cerr << "error: " << out.message << '\n';
    return 1;
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "--version") == 0) {
    std::cout << "flightsim " << flightsim::version_string() << '\n';
    return 0;
  }
  if (argc == 5 && std::strcmp(argv[1], "run") == 0 && std::strcmp(argv[3], "--out") == 0) {
    try {
      return cmd_run(argv[2], argv[4]);
    } catch (const std::exception& e) {
      std::cerr << "error: " << e.what() << '\n';
      return 2;
    }
  }
  return usage();
}
