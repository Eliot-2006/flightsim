#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <toml++/toml.hpp>

#include "flightsim/analysis/linearize.hpp"
#include "flightsim/analysis/modes.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/sim/aircraft_scenario.hpp"
#include "flightsim/sim/output.hpp"
#include "flightsim/sim/runner.hpp"
#include "flightsim/sim/scenario.hpp"
#include "flightsim/validation/convergence.hpp"
#include "flightsim/version.hpp"

namespace {

int usage() {
  std::cerr << "usage: flightsim --version\n"
               "       flightsim run <scenario.toml> --out <dir>\n"
               "       flightsim converge <study.toml> --out <dir>\n"
               "       flightsim trim <scenario.toml> --out <trim.json>\n"
               "       flightsim linearize <scenario.toml> --out <linear.json>\n";
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

// JSON strings use the existing TOML library's JSON formatter for escaping diagnostics and paths.
void json_string(std::ostream& out, const std::string& s) {
  out << toml::json_formatter{toml::value<std::string>{s}};
}
template <class Derived>
void json_vector(std::ostream& out, const Eigen::MatrixBase<Derived>& v) {
  out << '[';
  for (Eigen::Index i = 0; i < v.size(); ++i) {
    if (i) {
      out << ',';
    }
    out << v(i);
  }
  out << ']';
}
template <class Derived>
void json_matrix(std::ostream& out, const Eigen::MatrixBase<Derived>& m) {
  out << '[';
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    if (i) {
      out << ',';
    }
    json_vector(out, m.row(i));
  }
  out << ']';
}
int cmd_aircraft(const std::string& command, const std::string& scenario, const std::filesystem::path& dest) {
  using namespace flightsim;
  auto s = sim::load_aircraft_scenario(scenario);
  LookupMonitor trim_monitor;
  s.model.monitor = &trim_monitor;
  const auto tr = trim_level(s.model, s.trim);
  std::ostringstream out;
  out << std::setprecision(17) << std::boolalpha;
  out << "{\"xcg\":" << s.model.xcg << ",\"trim\":{\"ok\":" << tr.ok << ",\"alpha_rad\":" << tr.alpha
      << ",\"alpha_deg\":" << tr.alpha / units::deg << ",\"throttle\":" << tr.throttle
      << ",\"elevator_rad\":" << tr.de << ",\"elevator_deg\":" << tr.de / units::deg << ",\"residual\":";
  if (std::isfinite(tr.residual)) {
    out << tr.residual;
  } else {
    out << "null";
  }
  out << ",\"iterations\":" << tr.iterations << ",\"message\":";
  json_string(out, tr.message);
  out << ",\"lookup_out_of_domain\":" << trim_monitor.out_of_domain << ",\"first_lookup_violation\":";
  json_string(out, trim_monitor.first_violation);
  out << ",\"x17\":";
  json_vector(out, tr.x);
  out << '}';
  if (command == "linearize" && tr.ok) {
    const auto lm = linearize(s.model, tr);
    out << ",\"A\":";
    json_matrix(out, lm.A);
    out << ",\"B\":";
    json_matrix(out, lm.B);
    out << ",\"x0\":";
    json_vector(out, lm.x0);
    out << ",\"u0\":";
    json_vector(out, lm.u0);
    Eigen::Matrix4d lon_a, lat_a;
    Eigen::Matrix<double, 4, 2> lon_b;
    const int lon_ix[]{0, 1, 4, 7}, lat_ix[]{2, 3, 6, 8};
    for (int i = 0; i < 4; ++i) {
      const double row = i == 0 ? 1 / units::ft : 1;
      for (int j = 0; j < 4; ++j) {
        lon_a(i, j) = row * lm.A(lon_ix[i], lon_ix[j]) * (j == 0 ? units::ft : 1);
        lat_a(i, j) = lm.A(lat_ix[i], lat_ix[j]);
      }
      // Tabulated-engine input response: Pdot=0. State eigenvalues still hold P fixed.
      lon_b(i, 0) = row * (lm.B(lon_ix[i], 0) - lm.A(lon_ix[i], 12) * lm.B(12, 0) / lm.A(12, 12));
      lon_b(i, 1) = row * lm.B(lon_ix[i], 1) * units::deg;
    }
    out << ",\"sl_comparison\":{\"longitudinal_A\":";
    json_matrix(out, lon_a);
    out << ",\"longitudinal_B\":";
    json_matrix(out, lon_b);
    out << ",\"lateral_A\":";
    json_matrix(out, lat_a);
    out << '}';
    out << R"(,"state_names":["V","alpha","beta","phi","theta","psi","p","q","r","pN","pE","h","P"],
"state_units":["m/s","rad","rad","rad","rad","rad","rad/s","rad/s","rad/s","m","m","m","percent"],
"input_names":["throttle","elevator","aileron","rudder"],"input_units":["1","rad","rad","rad"],
"longitudinal_power":"held at trim; power state removed for modes","modes":[)";
    bool comma = false;
    for (const auto& mode : identify_modes(lm)) {
      if (comma) {
        out << ',';
      }
      comma = true;
      out << "{\"name\":";
      json_string(out, mode.name);
      out << ",\"real\":" << mode.eig.real() << ",\"imag\":" << mode.eig.imag() << ",\"wn\":" << mode.wn
          << ",\"zeta\":" << mode.zeta << ",\"tau\":" << mode.tau << '}';
    }
    out << "],\"chosen_epsilon\":6e-6,\"flat_region\":[";
    comma = false;
    for (double eps : {6e-5, 6e-7}) {
      const auto other = linearize(s.model, tr, eps);
      const double change = std::max((lm.A - other.A).cwiseAbs().maxCoeff() / lm.A.cwiseAbs().maxCoeff(),
                                     (lm.B - other.B).cwiseAbs().maxCoeff() / lm.B.cwiseAbs().maxCoeff());
      if (comma) {
        out << ',';
      }
      comma = true;
      out << "{\"epsilon\":" << eps << ",\"change\":" << change << '}';
    }
    out << "],\"sweep\":[";
    comma = false;
    for (const auto& point : linearization_sweep(s.model, tr)) {
      if (comma) {
        out << ',';
      }
      comma = true;
      out << "{\"epsilon\":" << point.eps << ",\"change\":" << point.change << '}';
    }
    // LN2 uses an interior altitude so the initial elevator lift transient stays in the thrust domain.
    auto response_spec = s.trim;
    response_spec.h = std::max(response_spec.h, 100 * units::ft);
    const auto response_trim = trim_level(s.model, response_spec);
    if (!response_trim.ok) {
      throw std::runtime_error("LN2 trim: " + response_trim.message);
    }
    const auto response_lm = linearize(s.model, response_trim);
    out << "],\"response_operating_point\":{\"h_m\":" << response_spec.h << ",\"V_mps\":" << response_spec.V
        << ",\"xcg\":" << response_spec.xcg << ",\"gamma_rad\":" << response_spec.gamma
        << "},\"response_x0\":";
    json_vector(out, response_lm.x0);
    out << ",\"responses\":[";
    comma = false;
    for (double step : {-0.5 * units::deg, -0.05 * units::deg}) {
      const auto response = elevator_response(s.model, response_trim, response_lm, step);
      if (comma) {
        out << ',';
      }
      comma = true;
      out << "{\"step_deg\":" << step / units::deg << ",\"max_error\":" << response.max_error
          << ",\"cells_changed\":" << response.cells_changed
          << ",\"out_of_domain\":" << response.out_of_domain
          << ",\"first_cell_change_time\":" << response.first_cell_change_time << ",\"first_violation\":";
      json_string(out, response.first_violation);
      out << ",\"max_abs_aileron\":" << response.max_abs_aileron
          << ",\"max_abs_rudder\":" << response.max_abs_rudder << ",\"cell_changes\":[";
      bool change_comma = false;
      for (const auto& change : response.cell_changes) {
        if (change_comma) {
          out << ',';
        }
        change_comma = true;
        out << "{\"table\":";
        json_string(out, change.table);
        out << ",\"from\":" << change.from << ",\"to\":" << change.to << ",\"t\":" << change.t
            << ",\"ignored\":" << change.ignored << ",\"reason\":";
        json_string(out, change.reason);
        out << '}';
      }
      out << "],\"points\":[";
      bool point_comma = false;
      for (const auto& point : response.points) {
        if (point_comma) {
          out << ',';
        }
        point_comma = true;
        out << "{\"t\":" << point.t << ",\"nonlinear\":";
        json_vector(out, point.nonlinear);
        out << ",\"linear\":";
        json_vector(out, point.linear);
        out << '}';
      }
      out << "]}";
    }
    out << ']';
  }
  out << "}\n";
  if (!dest.parent_path().empty()) {
    std::filesystem::create_directories(dest.parent_path());
  }
  std::ofstream file(dest);
  if (!file || !(file << out.str())) {
    throw std::runtime_error("cannot write " + dest.string());
  }
  std::cout << "wrote " << dest << " (" << tr.message << ")\n";
  return tr.ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "--version") == 0) {
    std::cout << "flightsim " << flightsim::version_string() << '\n';
    return 0;
  }
  if (argc == 5 &&
      (std::strcmp(argv[1], "run") == 0 || std::strcmp(argv[1], "converge") == 0 ||
       std::strcmp(argv[1], "trim") == 0 || std::strcmp(argv[1], "linearize") == 0) &&
      std::strcmp(argv[3], "--out") == 0) {
    try {
      if (std::strcmp(argv[1], "run") == 0) return cmd_run(argv[2], argv[4]);
      if (std::strcmp(argv[1], "trim") == 0 || std::strcmp(argv[1], "linearize") == 0) {
        return cmd_aircraft(argv[1], argv[2], argv[4]);
      }
      flightsim::validation::write_convergence_study(argv[2], argv[4]);
      std::cout << "wrote convergence.csv, orders.json and norm_conservation.csv to " << argv[4] << '\n';
      return 0;
    } catch (const std::exception& e) {
      std::cerr << "error: " << e.what() << '\n';
      return 2;
    }
  }
  return usage();
}
