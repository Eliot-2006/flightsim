#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "flightsim/aircraft/f16_data.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;
using units::deg;
using units::ft;
using units::lbf;
using units::slug;

namespace {

const std::filesystem::path kData = std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml";

const F16Data& real() {
  static const F16Data d = load_f16_data(kData);
  return d;
}

std::string slurp() {
  std::ifstream f(kData);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

std::filesystem::path write_tmp(const std::string& name, const std::string& text) {
  const auto p = std::filesystem::path(::testing::TempDir()) / name;
  std::ofstream(p) << text;
  return p;
}

// Both arguments are produced by the same IEEE operations only if the loader multiplies exactly as written,
// so these are tight: 1e-12 relative is ~4 orders above double eps and far below any unit-conversion mistake.
constexpr double kRel = 1e-12;

}  // namespace

TEST(F16Data, LoadsRealDataFile) {
  const auto& d = real();
  EXPECT_EQ(d.cx.x().size(), 12u);
  EXPECT_EQ(d.cx.y().size(), 5u);
  EXPECT_EQ(d.cm.x().size(), 12u);
  EXPECT_EQ(d.cm.y().size(), 5u);
  for (const Table2D* t : {&d.cl, &d.cn, &d.dlda, &d.dldr, &d.dnda, &d.dndr}) {
    EXPECT_EQ(t->x().size(), 12u) << t->name();
    EXPECT_EQ(t->y().size(), 7u) << t->name();
  }
  EXPECT_EQ(d.cz.x().size(), 12u);
  for (const Table1D* t : {&d.cxq, &d.cyr, &d.cyp, &d.czq, &d.clr, &d.clp, &d.cmq, &d.cnr, &d.cnp})
    EXPECT_EQ(t->x().size(), 12u) << t->name();  // 9 damping tables
  for (const Table2D* t : {&d.thrust_idle, &d.thrust_mil, &d.thrust_max}) {
    EXPECT_EQ(t->x().size(), 6u) << t->name();
    EXPECT_EQ(t->y().size(), 6u) << t->name();
  }
  // Coefficient values are untouched (TOML cx row alpha=0, elevator=0; resolved Cl; resolved CXq/CYr/Cz).
  EXPECT_EQ(d.cx(0.0, 0.0), -0.021);
  EXPECT_EQ(d.cl(15 * deg, 5 * deg), -0.019);
  EXPECT_EQ(d.cz(5 * deg), -0.416);
  EXPECT_EQ(d.cz(20 * deg), -1.366);
  EXPECT_EQ(d.cxq(-5 * deg), 0.11);
  EXPECT_EQ(d.cyr(40 * deg), -0.493);
  EXPECT_EQ(d.cmq(-5 * deg), -0.54);
  EXPECT_EQ(d.czq(0.0), -28.9);
  EXPECT_EQ(d.cnp(45 * deg), 0.15);
}

TEST(F16Data, MassPropsConvertedToSI) {
  const auto& d = real();
  const double mass = (1.0 / 0.00157) * slug;  // S&L uses reciprocal mass, not W/g
  EXPECT_NEAR(d.mass, mass, 1e-9 * mass);      // 1e-9 relative: required by the task spec
  const double i2 = slug * ft * ft;
  EXPECT_NEAR(d.Ixx, 9496 * i2, 1e-12 * 9496 * i2);
  EXPECT_NEAR(d.Iyy, 55814 * i2, kRel * 55814 * i2);
  EXPECT_NEAR(d.Izz, 63100 * i2, kRel * 63100 * i2);
  EXPECT_NEAR(d.Ixz, 982 * i2, kRel * 982 * i2);
  EXPECT_GT(d.Ixz, 0.0);  // positive product, not the tensor element
  EXPECT_NEAR(d.S, 300 * ft * ft, kRel * 300 * ft * ft);
  EXPECT_NEAR(d.b, 30 * ft, kRel * 30 * ft);
  EXPECT_NEAR(d.cbar, 11.32 * ft, kRel * 11.32 * ft);
  EXPECT_NEAR(d.hx, 160 * i2, kRel * 160 * i2);
  EXPECT_NEAR(d.weight_nominal, 20500 * lbf, kRel * 20500 * lbf);
  EXPECT_EQ(d.xcg_ref, 0.35);
  EXPECT_NEAR(d.g_sl, 32.17 * ft, kRel * 32.17 * ft);
  EXPECT_EQ(d.rtod, 57.29578);
}

TEST(F16Data, AngleBreakpointsInRadians) {
  const auto& d = real();
  EXPECT_NEAR(d.cx.x().front(), -10 * deg, 1e-15);  // 1e-15 abs: sub-ulp for O(0.1) rad
  EXPECT_NEAR(d.cx.x().back(), 45 * deg, 1e-15);
  EXPECT_NEAR(d.cx.y().front(), -24 * deg, 1e-15);
  EXPECT_NEAR(d.cx.y().back(), 24 * deg, 1e-15);
  EXPECT_NEAR(d.cl.y().front(), 0.0, 0.0);
  EXPECT_NEAR(d.cl.y().back(), 30 * deg, 1e-15);
  EXPECT_NEAR(d.dlda.y().front(), -30 * deg, 1e-15);  // effectiveness tables use signed beta
  EXPECT_NEAR(d.dndr.y().back(), 30 * deg, 1e-15);
  EXPECT_NEAR(d.cz.x().back(), 45 * deg, 1e-15);
  EXPECT_NEAR(d.cmq.x().front(), -10 * deg, 1e-15);
  EXPECT_NEAR(d.cz_beta_scale, 57.3 * deg, 1e-15);
  EXPECT_NEAR(d.cz_elevator_scale, 25 * deg, 1e-15);
  EXPECT_NEAR(d.aileron_norm, 20 * deg, 1e-15);
  EXPECT_NEAR(d.rudder_norm, 30 * deg, 1e-15);
  EXPECT_NEAR(d.cy_beta_gain, -0.02 / deg, 1e-12);  // per rad
  EXPECT_EQ(d.cy_aileron_gain, 0.021);
  EXPECT_EQ(d.cy_rudder_gain, 0.086);
  EXPECT_EQ(d.cz_elevator_gain, -0.19);
}

TEST(F16Data, ThrustInNewtons) {
  const auto& d = real();
  EXPECT_NEAR(d.thrust_idle.x().back(), 50000 * ft, 1e-9);  // altitude breakpoints in m
  EXPECT_EQ(d.thrust_idle.y().back(), 1.0);                 // Mach unchanged
  EXPECT_NEAR(d.thrust_idle(0.0, 0.0), 1060 * lbf, kRel * 1060 * lbf);
  EXPECT_NEAR(d.thrust_idle(10000 * ft, 0.6), -710 * lbf,
              kRel * 710 * lbf);  // resolved against NASA Table VI
  EXPECT_NEAR(d.thrust_mil(0.0, 0.0), 12680 * lbf, kRel * 12680 * lbf);
  EXPECT_GT(d.thrust_max(0.0, 0.0), d.thrust_mil(0.0, 0.0));
}

TEST(F16Data, EngineConstantsUnchangedAndActuatorsSI) {
  const auto& d = real();
  EXPECT_EQ(d.gear_breakpoint, 0.77);
  EXPECT_EQ(d.gear_low_slope, 64.94);
  EXPECT_EQ(d.gear_high_slope, 217.38);
  EXPECT_EQ(d.gear_high_intercept, -117.38);
  EXPECT_EQ(d.rtau_dp_low, 25.0);
  EXPECT_EQ(d.rtau_dp_high, 50.0);
  EXPECT_EQ(d.rtau_low, 1.0);
  EXPECT_EQ(d.rtau_high, 0.1);
  EXPECT_EQ(d.rtau_mid_intercept, 1.9);
  EXPECT_EQ(d.rtau_mid_slope, -0.036);
  EXPECT_EQ(d.power_threshold, 50.0);
  EXPECT_EQ(d.power_up_target, 60.0);
  EXPECT_EQ(d.power_down_target, 40.0);
  EXPECT_EQ(d.power_high_rate, 5.0);
  EXPECT_EQ(d.power_levels[0], 0.0);
  EXPECT_EQ(d.power_levels[1], 50.0);
  EXPECT_EQ(d.power_levels[2], 100.0);

  EXPECT_NEAR(d.elevator.pos_max, 25 * deg, 1e-15);
  EXPECT_NEAR(d.elevator.pos_min, -25 * deg, 1e-15);
  EXPECT_NEAR(d.aileron.pos_max, 21.5 * deg, 1e-15);
  EXPECT_NEAR(d.rudder.pos_max, 30 * deg, 1e-15);
  EXPECT_NEAR(d.elevator.rate_limit, 60 * deg, 1e-15);
  EXPECT_NEAR(d.aileron.rate_limit, 80 * deg, 1e-15);
  EXPECT_NEAR(d.rudder.rate_limit, 120 * deg, 1e-15);
  EXPECT_EQ(d.elevator.lag_rate, 20.2);
  EXPECT_EQ(d.throttle_min, 0.0);
  EXPECT_EQ(d.throttle_max, 1.0);
}

TEST(F16Data, SlAdcConstantsUnchangedImperial) {
  const auto& a = real().adc;
  EXPECT_EQ(a.rho0, 0.002377);
  EXPECT_EQ(a.temperature0, 519.0);
  EXPECT_EQ(a.temperature_high, 390.0);
  EXPECT_EQ(a.break_altitude, 35000.0);
  EXPECT_EQ(a.altitude_factor, 7.03e-06);
  EXPECT_EQ(a.density_exponent, 4.14);
  EXPECT_EQ(a.gamma, 1.4);
  EXPECT_EQ(a.gas_constant, 1716.3);
}

TEST(F16Data, RejectsMissingTable) {
  std::string s = slurp();
  const auto a = s.find("[aero.cn]");
  const auto b = s.find("\n[", a + 1);
  ASSERT_NE(a, std::string::npos);
  ASSERT_NE(b, std::string::npos);
  s.erase(a, b - a + 1);
  try {
    load_f16_data(write_tmp("no_cn.toml", s));
    FAIL() << "expected throw";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("aero.cn"), std::string::npos) << e.what();
  }
}

TEST(F16Data, RejectsUnknownTable) {
  const std::string s = slurp() + "\n[aero.bogus]\nvalue = 1\nsource = \"test\"\n";
  try {
    load_f16_data(write_tmp("extra.toml", s));
    FAIL() << "expected throw";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("aero.bogus"), std::string::npos) << e.what();
  }
}

TEST(F16Data, RejectsUnexpectedUnits) {
  std::string s = slurp();
  const auto a = s.find("[engine.idle]");
  const auto u = s.find("units = \"lbf\"", a);
  ASSERT_NE(u, std::string::npos);
  s.replace(u, 13, "units = \"N\"");
  try {
    load_f16_data(write_tmp("badunits.toml", s));
    FAIL() << "expected throw";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("engine.idle"), std::string::npos) << e.what();
  }
}

TEST(F16Data, RejectsMissingFile) {
  EXPECT_THROW(load_f16_data(write_tmp("absent_dir", "").parent_path() / "nope" / "x.toml"),
               std::runtime_error);
}
