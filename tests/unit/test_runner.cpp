#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "flightsim/sim/output.hpp"
#include "flightsim/sim/runner.hpp"
#include "flightsim/sim/scenario.hpp"
#include "flightsim/validation/analytic.hpp"

using namespace flightsim;
using namespace flightsim::sim;

namespace {

std::string scenario_text(const std::string& sim_extra, const std::string& dt, const std::string& tf,
                          const std::string& pos) {
  return "[sim]\nintegrator = \"rk4\"\ndt_s = " + dt + "\nt_final_s = " + tf + "\n" + sim_extra +
         "\n[body]\nmass_kg = 1.0\ninertia_kgm2 = [[1,0,0],[0,2,0],[0,0,3]]\n[initial]\nposition_m = " + pos +
         "\nvelocity_mps = [0,0,0]\neuler_deg = [0,0,0]\nomega_radps = [0,0,0]\n";
}

RigidBodyScenario load_text(const std::string& name, const std::string& text) {
  const auto p = std::filesystem::path(::testing::TempDir()) / name;
  std::ofstream(p) << text;
  return load_rigid_body_scenario(p);
}

std::string slurp(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

RigidBodyScenario linear_drag_example() {
  return load_rigid_body_scenario(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "scenarios" /
                                  "rb_linear_drag.toml");
}

}  // namespace

TEST(Runner, NonFiniteStops) {
  auto s = load_text("nan.toml",
                     scenario_text("", "0.125", "1.0", "[0,0,0]"));  // dt exact in binary: boundaries k/8
  s.forces.push_back([](double t, const Vec<13>&) {
    ForceMoment fm;
    if (t >= 0.5) fm.F_N = Vec3::Constant(std::nan(""));
    return fm;
  });
  const RunOutput out = run(s);
  EXPECT_EQ(out.reason, StopReason::NonFinite);
  // RK4's last stage of the step 0.375 -> 0.5 samples the force at t = 0.5, so the first non-finite boundary
  // is 0.5.
  EXPECT_EQ(out.t_stop, 0.5);
  EXPECT_NE(out.message.find("non-finite"), std::string::npos) << out.message;
  EXPECT_NE(out.message.find("t=0.5"), std::string::npos) << out.message;
  EXPECT_NE(out.message.find("nan"), std::string::npos) << out.message;  // lists the offending state
  ASSERT_FALSE(out.samples.empty());
  EXPECT_LT(out.samples.back().t, 0.5);  // only finite samples are recorded
  for (const auto& smp : out.samples) EXPECT_TRUE(smp.x.allFinite());
}

TEST(Runner, GroundContactStopsWhenHBelowZero) {
  // Drop from h = 10 m: h = 10 - g t^2/2 reaches 0 at t = sqrt(20/g) = 1.4279 s -> first boundary at 1.43 s.
  const auto s = load_text("ground.toml", scenario_text("stop_at_ground = true", "0.01", "5.0", "[0,0,-10]"));
  const RunOutput out = run(s);
  EXPECT_EQ(out.reason, StopReason::GroundContact);
  EXPECT_NEAR(out.t_stop, 1.43, 1e-9);  // 1e-9: t_k = k*dt carries ~1e-16 roundoff
  ASSERT_GE(out.samples.size(), 2u);
  EXPECT_LE(-out.samples.back().x(2), 0.0);
  EXPECT_GT(-out.samples[out.samples.size() - 2].x(2), 0.0);
  EXPECT_EQ(out.samples.back().t, out.t_stop);

  // Without the flag the same run goes through the ground to the final time.
  const auto s2 = load_text("noground.toml", scenario_text("", "0.01", "5.0", "[0,0,-10]"));
  const RunOutput out2 = run(s2);
  EXPECT_EQ(out2.reason, StopReason::FinalTime);
  EXPECT_EQ(out2.t_stop, 5.0);
}

TEST(Runner, DeterministicBitwise) {
  const auto s = linear_drag_example();
  const auto dir = std::filesystem::path(::testing::TempDir());
  const RunOutput a = run(s), b = run(s);
  write_csv(dir / "det_a.csv", s, a);
  write_csv(dir / "det_b.csv", s, b);
  const std::string ca = slurp(dir / "det_a.csv"), cb = slurp(dir / "det_b.csv");
  EXPECT_FALSE(ca.empty());
  EXPECT_EQ(ca, cb);  // same binary + scenario => identical bytes (PRD 9.6)
}

TEST(Runner, SamplingAndAgreementWithAnalytic) {
  const auto s = linear_drag_example();
  const RunOutput out = run(s);
  EXPECT_EQ(out.reason, StopReason::FinalTime);
  EXPECT_EQ(out.t_stop, s.run.t_final);
  const long n = std::lround(s.run.t_final / s.run.dt);
  EXPECT_EQ(static_cast<long>(out.samples.size()), n / s.run.output_every + 1);
  EXPECT_EQ(out.nfev, 4 * n);  // RK4: 4 derivative evaluations per step
  EXPECT_EQ(out.samples.front().t, 0.0);
  EXPECT_EQ(out.samples.back().t, s.run.t_final);

  // Wiring check against the T2 closed form (same parameters as scenarios/rb_linear_drag.toml).
  const auto ex =
      validation::exact_linear_drag(2.0, 0.5, s.g_N, s.x0.head<3>(), s.x0.segment<3>(3), s.run.t_final);
  // RK4 at dt = 0.01 on this smooth problem: error ~ 1e-13; 1e-9 leaves four orders of margin.
  EXPECT_LE((out.samples.back().x.head<3>() - ex.p).norm(), 1e-9);
  EXPECT_LE((out.samples.back().x.segment<3>(3) - ex.v).norm(), 1e-9);
}

TEST(Runner, WritesCsvAndMetadata) {
  const auto s = linear_drag_example();
  const RunOutput out = run(s);
  const auto dir = std::filesystem::path(::testing::TempDir());
  write_csv(dir / "traj.csv", s, out);
  write_metadata(dir / "meta.json", s, out);
  const std::string csv = slurp(dir / "traj.csv");
  EXPECT_EQ(csv.substr(0, 10), "t_s,pN_m,p");
  EXPECT_EQ(std::count(csv.begin(), csv.end(), '\n'), static_cast<long>(out.samples.size()) + 1);
  for (const char* col : {"phi_rad", "psi_rad", "h_m", "V_mps", "alpha_rad", "E_J", "HN_x_kgm2ps", "qnorm"})
    EXPECT_NE(csv.find(col), std::string::npos) << col;
  const std::string meta = slurp(dir / "meta.json");
  for (const char* key :
       {"git_sha", "build_type", "compiler", "scenario_hash", "integrator", "dt_s", "stop_reason", "nfev"})
    EXPECT_NE(meta.find(key), std::string::npos) << key;
  EXPECT_NE(meta.find("rk4"), std::string::npos);
}
