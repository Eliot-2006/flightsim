#include <gtest/gtest.h>

#include <toml++/toml.hpp>

#include "flightsim/analysis/trim.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/sim/aircraft_scenario.hpp"

using namespace flightsim;
namespace {
const std::filesystem::path root = FLIGHTSIM_SOURCE_DIR;
}

TEST(Trim, TrimConvergesOnSyntheticQuadratic) {
  const auto s = solve_trim_system(
      [](const Vec3& z) -> Vec3 { return Vec3(z(0) * z(0) - 4, z(1) * z(1) - 9, z(2) * z(2) - 16); },
      Vec3(1, 2, 3));
  ASSERT_TRUE(s.ok);
  // Task 5.1 residual convergence threshold; positive-root error is smaller than residual.
  EXPECT_LT((s.z - Vec3(2, 3, 4)).lpNorm<Eigen::Infinity>(), 1e-10);
  EXPECT_LT(s.residual, 1e-10);
}

TEST(Trim, TrimReportsFailureOutsideEnvelope) {
  auto scenario = sim::load_aircraft_scenario(root / "scenarios/f16_trim_502_xcg035.toml");
  scenario.trim.V = 60 * units::ft;
  const auto tr = trim_level(scenario.model, scenario.trim);
  EXPECT_FALSE(tr.ok);
  EXPECT_FALSE(tr.message.empty());
}

TEST(AircraftScenario, UnitsPathsAndStrictSchema) {
  auto t = toml::parse(
      "[aircraft]\ndata='data/f16/f16_sl.toml'\nxcg=0.3\natmosphere='sl_adc'\n"
      "gravity='sl'\nrotor=true\nactuators=false\n[trim]\nV_fps=502\nh_ft=0\ngamma_deg=0\n");
  auto m = sim::load_aircraft(t, root);
  EXPECT_DOUBLE_EQ(m.grav.g0, m.data.g_sl);
  EXPECT_TRUE(m.include_rotor);
  EXPECT_FALSE(m.actuators);
  const auto spec = sim::load_trim_spec(t, m.xcg);
  EXPECT_DOUBLE_EQ(spec.V, 502 * units::ft);
  t["trim"].as_table()->insert("V_mps", 153.0);
  EXPECT_THROW(sim::load_trim_spec(t, m.xcg), std::runtime_error);
  t["aircraft"].as_table()->insert("typo", true);
  EXPECT_THROW(sim::load_aircraft(t, root), std::runtime_error);
}

TEST(AircraftScenario, GravityAndUnknownKeys) {
  const std::string prefix =
      "[aircraft]\ndata='data/f16/f16_sl.toml'\nxcg=0.3\natmosphere='sl_adc'\n"
      "rotor=true\nactuators=false\n";
  auto t = toml::parse(prefix + "gravity={g0=9.8,R=inf}\n[trim]\nV_mps=153\nh_m=0\ngamma_deg=0\n");
  const auto m = sim::load_aircraft(t, root);
  EXPECT_DOUBLE_EQ(m.grav.g0, 9.8);
  EXPECT_TRUE(std::isinf(m.grav.R));
  EXPECT_DOUBLE_EQ(sim::load_trim_spec(t, .3).V, 153);
  for (const auto& extra : {"\nunknown=1\n", "\n[extra]\nunknown=1\n"}) {
    auto bad = toml::parse(prefix + "gravity='sl'\n[trim]\nV_mps=153\nh_m=0\ngamma_deg=0\n" + extra);
    EXPECT_THROW(sim::load_trim_spec(bad, .3), std::runtime_error);
  }
  auto bad = toml::parse(prefix + "gravity={g0=9.8,R=-1}\n");
  EXPECT_THROW(sim::load_aircraft(bad, root), std::runtime_error);
}

TEST(Trim, TrimRejectsOutOfTableEngine) {
  auto s = sim::load_aircraft_scenario(root / "scenarios/f16_trim_502_xcg035.toml");
  for (const auto& spec :
       {TrimSpec{502 * units::ft, -100 * units::ft, 0, .35}, TrimSpec{1150 * units::ft, 0, 0, .35}}) {
    const auto tr = trim_level(s.model, spec);
    EXPECT_FALSE(tr.ok);
    EXPECT_NE(tr.message.find("outside"), std::string::npos);
  }
}

TEST(Trim, ActuatorEnabledTrimIsEquilibriumForActuatorAndPowerStates) {
  // The control scenario enables actuators; trim solves on ideal surfaces, then sets the actuator states
  // to the trimmed positions (zero actuator derivative).
  auto s = sim::load_aircraft_scenario(root / "scenarios/control/f16_ap_trim.toml");
  ASSERT_TRUE(s.model.actuators);
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok) << tr.message;
  EXPECT_DOUBLE_EQ(tr.x(ac::kDe), tr.de);
  EXPECT_DOUBLE_EQ(tr.x(ac::kDa), 0);
  EXPECT_DOUBLE_EQ(tr.x(ac::kDr), 0);
  const Vec<17> xd = s.model.deriv(0, tr.x, tr.u);
  // Actuator and power derivatives vanish to roundoff; 1e-12 (1/s-scale units, state magnitudes <= 1 rad or
  // ~11 %) is ~4 orders above eps and ~2 below the 1e-10 trim residual that sets the rigid-body entries.
  for (int i : {ac::kPower, ac::kDe, ac::kDa, ac::kDr}) {
    EXPECT_NEAR(xd(i), 0, 1e-12) << "state " << i;
  }
  // Same trim as the ideal-surface model: the actuator-enabled model changes nothing but the extra states.
  auto ideal = s;
  ideal.model.actuators = false;
  const auto tr_ideal = trim_level(ideal.model, ideal.trim);
  EXPECT_EQ(tr.alpha, tr_ideal.alpha);
  EXPECT_EQ(tr.throttle, tr_ideal.throttle);
  EXPECT_EQ(tr.de, tr_ideal.de);
}
