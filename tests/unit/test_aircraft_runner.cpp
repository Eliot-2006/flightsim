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
  { std::ofstream(path) << text; }
  return sim::load_aircraft_scenario(path);
}
}  // namespace

TEST(AircraftRunner, RejectsUnalignedInputEvent) {
  // dt = 0.01 s: an event at 0.505 s falls between steps.
  EXPECT_THROW(
      load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0.505\nt_end_s=1.0\nvalue=1\n")),
      std::runtime_error);
  EXPECT_NO_THROW(
      load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0.5\nt_end_s=1.0\nvalue=1\n")));
}

TEST(AircraftRunner, DoubletAppliedOnSchedule) {
  const auto s = load(scenario_text("3\n",
                                    "[[inputs]]\nchannel='de'\nt_start_s=1\nt_end_s=2\nvalue=1\n"
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
  const auto out = sim::run_aircraft(
      load(scenario_text("2\n", "[[inputs]]\nchannel='de'\nt_start_s=0\nt_end_s=2\nvalue=60\n")));
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

namespace {
// Closed-loop scenario text: actuators on, 50 Hz autopilot with the committed gains, then `extra`.
std::string ap_text(const std::string& sim_extra, const std::string& extra, double rate_hz = 50,
                    double dt = 0.01) {
  return "[aircraft]\ndata='" + root +
         "/data/f16/f16_sl.toml'\nxcg=0.35\natmosphere='sl_adc'\ngravity='sl'\nrotor=true\nactuators=true\n"
         "[trim]\nV_fps=502\nh_ft=10000\ngamma_deg=0\n[sim]\nintegrator='rk4'\ndt_s=" +
         std::to_string(dt) + "\nt_final_s=" + sim_extra +
         "\n[autopilot]\nrate_hz=" + std::to_string(rate_hz) + "\ngains='" + root +
         "/scenarios/control/autopilot_gains.toml'\nantiwindup=true\n" + extra;
}
const char* kAltCmd = "[[commands]]\nkind='altitude'\nt_s=1\nvalue_ft=10100\n";
}  // namespace

TEST(AircraftRunnerAutopilot, RejectsDtThatDoesNotDivideTheControllerPeriod) {
  // 1/50 s is not a multiple of 0.03 s.
  try {
    load(ap_text("3.0\n", "", 50, 0.03));
    FAIL() << "expected a load error";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("dt_s"), std::string::npos) << e.what();
  }
  EXPECT_NO_THROW(load(ap_text("3.0\n", "")));            // 0.01 divides 0.02
  EXPECT_NO_THROW(load(ap_text("3.0\n", "", 50, 0.02)));  // one step per sample
}

TEST(AircraftRunnerAutopilot, RejectsRateThatDiffersFromTheDesignRate) {
  // The gains were designed (Pade delay, filter, tracking times) for 50 Hz.
  EXPECT_THROW(load(ap_text("3.0\n", "", 25)), std::runtime_error);
}

TEST(AircraftRunnerAutopilot, RejectsInconsistentTables) {
  EXPECT_THROW(load(ap_text("3.0\n", "[[inputs]]\nchannel='de'\nt_start_s=0\nt_end_s=1\nvalue=1\n")),
               std::runtime_error);  // open-loop events and a controller both drive the surfaces
  EXPECT_THROW(load(ap_text("3.0\n", "[[commands]]\nkind='heading'\nt_s=1\nvalue_ft=1\n")),
               std::runtime_error);
  EXPECT_THROW(load(ap_text("3.0\n", "[[commands]]\nkind='altitude'\nt_s=1.005\nvalue_ft=10100\n")),
               std::runtime_error);  // not on a step boundary
  EXPECT_THROW(load(ap_text("3.0\n", "[[commands]]\nkind='altitude'\nt_s=1\nvalue_ft=1\nvalue_m=1\n")),
               std::runtime_error);  // exactly one unit
  EXPECT_THROW(load(ap_text("3.0\n", "[[commands]]\nkind='altitude'\nt_s=1\nvalue_fps=3\n")),
               std::runtime_error);  // airspeed unit on an altitude command
  EXPECT_THROW(load(scenario_text("3.0\n", kAltCmd)), std::runtime_error);  // commands without [autopilot]
  EXPECT_THROW(load(ap_text("3.0\n", "[autopilot]\n")), std::runtime_error);
}

TEST(AircraftRunnerAutopilot, EngagedWithoutCommandsHoldsTrim) {
  const auto s = load(ap_text("5.0\n", ""));
  const auto out = sim::run_aircraft(s);
  ASSERT_EQ(out.reason, sim::StopReason::FinalTime);
  EXPECT_FALSE(out.out_of_domain) << out.first_violation;
  // Bumpless engage: the very first command is the trim command, exactly, and the actuator states start on
  // it.
  EXPECT_EQ(out.u.front().de, out.u_trim.de);
  EXPECT_EQ(out.u.front().throttle, out.u_trim.throttle);
  EXPECT_EQ(out.x.front()(ac::kDe), out.u_trim.de);
  const Vec<13> a = x13_from_x17(s.model, out.x.front()), b = x13_from_x17(s.model, out.x.back());
  // Trim residual < 1e-10 (scaled) over 5 s; the loops add feedback on errors of that size only.
  EXPECT_NEAR(b(11), a(11), 1e-4);  // altitude [m]
  EXPECT_NEAR(b(0), a(0), 1e-6);    // airspeed [m/s]
  EXPECT_NEAR(out.u.back().de, out.u_trim.de, 1e-6);
  EXPECT_NEAR(out.u.back().throttle, out.u_trim.throttle, 1e-6);
}

TEST(AircraftRunnerAutopilot, OutputsAreHeldBetweenControllerSamples) {
  const auto out = sim::run_aircraft(load(ap_text("3.0\n", kAltCmd)));
  ASSERT_EQ(out.t.size(), 301u);
  bool changed = false;
  for (std::size_t i = 1; i < out.u.size(); ++i) {
    // dt = 0.01 s, controller 50 Hz: a new command can only appear at even steps (samples 0, 2, 4, ...).
    if (i % 2 == 1) {
      EXPECT_EQ(out.u[i].de, out.u[i - 1].de) << i;
      EXPECT_EQ(out.u[i].throttle, out.u[i - 1].throttle) << i;
    } else if (out.u[i].de != out.u[i - 1].de) {
      changed = true;
    }
  }
  EXPECT_TRUE(changed);
  // Before the command (t < 1 s) the loop holds trim; the +100 ft command at t = 1 s moves the elevator
  // trailing-edge-up within the next samples.
  // 1e-9 rad: the trim residual (< 1e-10 scaled) is fed back by the loops (observed drift ~5e-13 rad).
  EXPECT_NEAR(out.u[99].de, out.u[0].de, 1e-9);
  EXPECT_LT(out.u[110].de, out.u[0].de);
}

TEST(AircraftRunnerAutopilot, TelemetryAndMetadataCarryTheCommand) {
  const auto s = load(ap_text("3.0\n", kAltCmd));
  const auto out = sim::run_aircraft(s);
  ASSERT_EQ(out.ap.size(), out.t.size());
  const double h_trim = x13_from_x17(s.model, out.x.front())(11);
  EXPECT_DOUBLE_EQ(out.ap[50].h_cmd, h_trim);
  EXPECT_NEAR(out.ap[150].h_cmd, 10100 * units::ft, 1e-9);
  const auto csv = std::filesystem::temp_directory_path() / "flightsim_test_ap.csv";
  sim::write_aircraft_csv(csv, s, out);
  std::ifstream f(csv);
  std::string header;
  std::getline(f, header);
  for (const char* col : {"h_cmd_m", "V_cmd_mps", "theta_cmd_rad", "q_cmd_radps"}) {
    EXPECT_NE(header.find(col), std::string::npos) << col;
  }
}

TEST(AircraftRunnerAutopilot, DeterministicBitwise) {
  const auto s = load(ap_text("3.0\n", kAltCmd));
  const auto a = sim::run_aircraft(s), b = sim::run_aircraft(s);
  ASSERT_EQ(a.x.size(), b.x.size());
  EXPECT_EQ(0, std::memcmp(a.x.back().data(), b.x.back().data(), sizeof(double) * 17));
}

TEST(AircraftRunnerWind, ConstantWindStepIsAppliedFromItsStepBoundary) {
  // Headwind 10 m/s (wind from the north, aircraft heading north) switched on at t = 1 s.
  const auto s = load(scenario_text("2.0\n", "[[wind.steps]]\nt_s=1\ndelta_mps=[-10.0, 0.0, 0.0]\n"));
  EXPECT_EQ(sim::constant_wind_at(s, 0.99)(0), 0.0);
  EXPECT_EQ(sim::constant_wind_at(s, 1.0)(0), -10.0);
  const auto out = sim::run_aircraft(s);
  const double before = x13_from_x17(s.model, out.x[99], 0.99)(0);
  auto windy = s.model;
  windy.wind.constant_N = Vec3(-10, 0, 0);
  const double after = x13_from_x17(windy, out.x[100], 1.0)(0);
  // Inertial velocity is continuous across the step, so the air-relative speed jumps by the wind change
  // (less the ~0.01 s of acceleration, <= 0.1 m/s).
  EXPECT_NEAR(after - before, 10.0, 0.1);
  EXPECT_THROW(load(scenario_text("2.0\n", "[[wind.steps]]\nt_s=1.005\ndelta_mps=[-10.0, 0.0, 0.0]\n")),
               std::runtime_error);
}
