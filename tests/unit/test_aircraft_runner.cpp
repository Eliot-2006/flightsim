#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <fstream>
#include <toml++/toml.hpp>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/sim/aircraft_runner.hpp"

using namespace flightsim;
namespace {
const std::string root = FLIGHTSIM_SOURCE_DIR;

// Scenario text with an absolute data path; `extra` is appended (inputs, wind, ...).
std::string scenario_text(const std::string& sim_extra, const std::string& extra,
                          const std::string& atmosphere = "sl_adc") {
  return "[aircraft]\ndata='" + root + "/data/f16/f16_sl.toml'\nxcg=0.35\natmosphere='" + atmosphere +
         "'\ngravity='sl'\nrotor=true\nactuators=false\n[trim]\nV_fps=502\nh_ft=10000\ngamma_deg=0\n"
         "[sim]\nintegrator='rk4'\ndt_s=0.01\nt_final_s=" +
         sim_extra + "\n" + extra;
}
sim::AircraftScenario load(const std::string& text) {
  const auto path = std::filesystem::temp_directory_path() / "flightsim_test_aircraft.toml";
  {
    std::ofstream(path) << text;
  }
  return sim::load_aircraft_scenario(path);
}
}  // namespace

TEST(AircraftRunner, RejectsUnalignedInputEvent) {
  // dt = 0.01 s: an event at 0.505 s falls between steps.
  EXPECT_THROW(load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0.505\nt_end_s=1.0\nvalue=1\n")),
               std::runtime_error);
  EXPECT_NO_THROW(load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0.5\nt_end_s=1.0\nvalue=1\n")));
}

TEST(AircraftRunner, DoubletAppliedOnSchedule) {
  const auto s = load(scenario_text(
      "3\n", "[[inputs]]\nchannel='de'\nt_start_s=1\nt_end_s=2\nvalue=1\n"
             "[[inputs]]\nchannel='de'\nt_start_s=2\nt_end_s=3\nvalue=-1\n"));
  const auto out = sim::run_aircraft(s);
  ASSERT_EQ(out.t.size(), 301u);
  for (std::size_t i = 0; i < out.t.size(); ++i) {
    const double t = out.t[i];
    const double expected = (t >= 1 && t < 2) ? 1.0 : (t >= 2 && t < 3) ? -1.0 : 0.0;
    // t is k*dt in double, so compare with half a step of slack on the schedule boundaries.
    if (std::abs(t - 1) < 0.005 || std::abs(t - 2) < 0.005 || std::abs(t - 3) < 0.005) {
      continue;
    }
    EXPECT_NEAR(out.u[i].de - out.u_trim.de, expected * units::deg, 1e-15) << "t=" << t;
  }
  // The boundary samples themselves: input at t=1 is +1 deg (event active from its first step).
  EXPECT_NEAR(out.u[100].de - out.u_trim.de, units::deg, 1e-15);
  EXPECT_NEAR(out.u[200].de - out.u_trim.de, -units::deg, 1e-15);
}

TEST(AircraftRunner, FromTrimStaysTrimmedWithNoInputs) {
  const auto out = sim::run_aircraft(load(scenario_text("2\n", "")));
  const auto model = load(scenario_text("2\n", "")).model;
  const double v0 = x13_from_x17(model, out.x.front())(0);
  const double v1 = x13_from_x17(model, out.x.back())(0);
  // Trim residual is < 1e-10 (scaled by g); 2 s of drift is ~1e-9 m/s, so 1e-6 m/s is a loose bound.
  EXPECT_LT(std::abs(v1 - v0), 1e-6);
  EXPECT_FALSE(out.out_of_domain);
}

TEST(AircraftRunner, OutOfDomainReported) {
  // +60 deg elevator step: alpha and the elevator table are driven out of their domains.
  const auto out =
      sim::run_aircraft(load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0\nt_end_s=2\nvalue=60\n")));
  EXPECT_TRUE(out.out_of_domain);
  EXPECT_FALSE(out.first_violation.empty());
}

TEST(AircraftRunner, DeterministicBitwise) {
  const auto s = load(scenario_text("2\n", "[[inputs]]\nchannel='da'\nt_start_s=0.5\nt_end_s=1\nvalue=2\n"));
  const auto a = sim::run_aircraft(s), b = sim::run_aircraft(s);
  ASSERT_EQ(a.x.size(), b.x.size());
  EXPECT_EQ(0, std::memcmp(a.x.back().data(), b.x.back().data(), sizeof(double) * 17));
}

TEST(AircraftRunner, CsvAndConstantWindParse) {
  auto s = load(scenario_text("1\n", "[wind]\nconstant_mps=[3.0, -2.0, 0.0]\n"));
  EXPECT_DOUBLE_EQ(s.model.wind.constant_N(1), -2.0);
  const auto out = sim::run_aircraft(s);
  const auto csv = std::filesystem::temp_directory_path() / "flightsim_test_aircraft.csv";
  sim::write_aircraft_csv(csv, s, out);
  std::ifstream f(csv);
  std::string header;
  std::getline(f, header);
  EXPECT_NE(header.find("V_fps"), std::string::npos);
  EXPECT_NE(header.find("alpha_rad"), std::string::npos);
}
