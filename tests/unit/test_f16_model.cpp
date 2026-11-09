#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <toml++/toml.hpp>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;

namespace {
F16Model model() {
  F16Model m{load_f16_data(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml")};
  m.grav.g0 = m.data.g_sl;
  return m;
}

Vec<13> flight_state() {
  Vec<13> x = Vec<13>::Zero();
  x(0) = 152.4;
  x(1) = x(4) = 5 * units::deg;
  x(11) = 3048;
  x(12) = 50;
  return x;
}

double value(const toml::node_view<const toml::node>& table, const char* key) {
  return table[key]["value"].value<double>().value();
}

void check_reference(const char* name, double actual, double expected) {
  const double error = std::abs((actual - expected) / expected);
  std::cout << std::setprecision(12) << "A1 " << name << " actual=" << actual << " expected=" << expected
            << " relative_error=" << error << '\n';
  // PRD 8.5 / Task 4.5: published S&L values have only 4-6 printed significant figures.
  EXPECT_LE(error, 1e-3) << name;
  // Investigated list (docs/results.md): p_dot/r_dot use rounded S&L c1..c9 in the publication;
  // production derives them from the full inertia. A1 below verifies that those rounded numbers close
  // the discrepancy. The binding 1e-3 check above applies to every component, including these two.
  const bool rounded_inertia = std::string_view(name) == "p_dot" || std::string_view(name) == "r_dot";
  if (error > 1e-4 && !rounded_inertia) {
    ADD_FAILURE() << name << ": investigate relative error " << error;
  }
}
}  // namespace

TEST(F16Model, A1CheckCase) {
  auto m = model();
  const auto refs = toml::parse_file(
      (std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/reference_values.toml").string());
  const auto a1 = refs["A1"];
  m.xcg = value(a1, "xcg");
  const std::array<const char*, 13> names = {"vt", "alpha", "beta", "phi", "theta", "psi",  "p",
                                             "q",  "r",     "pn",   "pe",  "h",     "power"};
  Vec<13> x;
  for (int i = 0; i < 13; ++i) {
    x(i) = value(a1, names[static_cast<std::size_t>(i)]);
    if (i == 0 || (i >= 9 && i <= 11)) x(i) *= units::ft;
  }
  const F16Input u{value(a1, "throttle"), value(a1, "elevator") * units::deg,
                   value(a1, "aileron") * units::deg, value(a1, "rudder") * units::deg};
  Vec<13> d = f13(m, x, u);
  for (int i = 0; i < 13; ++i) {
    if (i == 0 || (i >= 9 && i <= 11)) d(i) /= units::ft;
    const std::string name = std::string(names[static_cast<std::size_t>(i)]) + "_dot";
    check_reference(name.c_str(), d(i), value(a1["expected"], name.c_str()));
  }
  const auto x17 = x17_from_x13(m, x, u);
  const auto d17 = m.deriv(0, x17, u);
  const auto fm = m.evaluate(0, x17, u).fm;
  std::cout << "A1 moments_lbf_ft=" << (fm.M_B / (units::lbf * units::ft)).transpose() << '\n';
  const Vec3 moment = fm.M_B / (units::lbf * units::ft);
  const double p = x(6), q = x(7), r = x(8);
  // Published S&L rounded inertia combinations, diagnostic only (SOURCES.md "Mass and inertia").
  const Vec3 rounded(
      (0.02755 * p - 0.770 * r + 1.642e-6 * 160) * q + 1.055e-4 * moment.x() + 1.642e-6 * moment.z(),
      0.9604 * p * r - 0.01759 * (p * p - r * r) + 1.792e-5 * (moment.y() - 160 * r),
      (-0.7336 * p - 0.02755 * r + 1.587e-5 * 160) * q + 1.642e-6 * moment.x() + 1.587e-5 * moment.z());
  std::cout << "A1 rounded_c1_c9_rates=" << rounded.transpose() << '\n';
  for (int i = 0; i < 3; ++i) {
    const auto key = std::array<const char*, 3>{"p_dot", "q_dot", "r_dot"}[static_cast<std::size_t>(i)];
    const double expected = value(a1["expected"], key);
    // 1e-6 relative diagnostic: the rounded coefficients reproduce the printed derivatives (6-7 digits).
    EXPECT_NEAR(rounded(i), expected, 1e-6 * std::abs(expected)) << key;
  }
  const Mat3 C_BN = dcm_from_quat(rb::quat(x17)).transpose();
  const Vec3 v = C_BN * rb::vel(x17);
  const Vec3 acceleration = C_BN * rb::vel(d17) - rb::omega(x17).cross(v);
  for (int i = 0; i < 3; ++i) {
    const char* key = std::array<const char*, 3>{"u_dot", "v_dot", "w_dot"}[static_cast<std::size_t>(i)];
    check_reference(key, acceleration(i) / units::ft, value(a1["body_axes"], key));
  }
}

TEST(F16Model, PositiveElevatorNoseDown) {
  const auto m = model();
  const auto x = x17_from_x13(m, flight_state(), {});
  EXPECT_LT(m.evaluate(0, x, {0, units::deg, 0, 0}).fm.M_B.y(), m.evaluate(0, x, {}).fm.M_B.y());
}

TEST(F16Model, PositiveAileronLeftRoll) {
  const auto m = model();
  const auto x = x17_from_x13(m, flight_state(), {});
  EXPECT_LT(m.evaluate(0, x, {0, 0, units::deg, 0}).fm.M_B.x(), m.evaluate(0, x, {}).fm.M_B.x());
}

TEST(F16Model, PositiveRudderLeftYaw) {
  const auto m = model();
  const auto x = x17_from_x13(m, flight_state(), {});
  EXPECT_LT(m.evaluate(0, x, {0, 0, 0, units::deg}).fm.M_B.z(), m.evaluate(0, x, {}).fm.M_B.z());
}

TEST(F16Model, PitchDampingNegative) {
  const auto m = model();
  auto x = x17_from_x13(m, flight_state(), {});
  const double baseline = m.evaluate(0, x, {}).fm.M_B.y();
  x(rb::kOmega + 1) = 0.1;
  EXPECT_LT(m.evaluate(0, x, {}).fm.M_B.y(), baseline);
}

TEST(F16Model, ThrustIncreasesWithThrottle) {
  const auto m = model();
  double previous = -std::numeric_limits<double>::infinity();
  for (double throttle : {0.0, 0.3, 0.7, 0.9, 1.0}) {
    auto x = x17_from_x13(m, flight_state(), {});
    const F16Input u{throttle, 0, 0, 0};
    for (int i = 0; i < 4000; ++i) x(ac::kPower) += 0.01 * m.deriv(0, x, u)(ac::kPower);
    const auto e = m.evaluate(0, x, u);
    // 40 s is at least 40 final time constants (1/s); 1e-10 percent covers integration roundoff.
    EXPECT_NEAR(x(ac::kPower), e.power_cmd, 1e-10);
    EXPECT_GT(e.thrust, previous);
    previous = e.thrust;
  }
}

TEST(F16Model, EnginePowerLagTransitions) {
  const auto m = model();
  struct Case {
    double power, throttle, derivative;
  };
  const Case cases[] = {{10, 0.3, 9.482}, {40, 1, 20},       {49, 1, 11},   {50, 0, -50},
                        {70, 1, 150},     {70, 0, -150},     {0, 1, 6},     {0, 0.5, 23.7381676},
                        {0, 0.77, 6},     {90, 0.9, -58.69}, {40, 0, -18.4}};
  auto x = x17_from_x13(m, flight_state(), {});
  for (const auto& c : cases) {
    x(ac::kPower) = c.power;
    // Hand evaluation of S&L gearing/rtau/transition branches; 1e-10 percent/s for decimal rounding.
    EXPECT_NEAR(m.deriv(0, x, {c.throttle, 0, 0, 0})(ac::kPower), c.derivative, 1e-10)
        << "power=" << c.power << " throttle=" << c.throttle;
  }
  x(ac::kPower) = 0;
  EXPECT_DOUBLE_EQ(m.evaluate(0, x, {-1, 0, 0, 0}).power_cmd, 0);
  EXPECT_DOUBLE_EQ(m.evaluate(0, x, {2, 0, 0, 0}).power_cmd, 100);
}

TEST(F16Model, ActuatorRateAndPositionLimits) {
  auto m = model();
  m.actuators = true;
  const F16Input u{0.5, 1, -1, 1};
  auto x = x17_from_x13(m, flight_state(), {});
  const auto d = m.deriv(0, x, u);
  EXPECT_DOUBLE_EQ(d(ac::kDe), m.data.elevator.rate_limit);
  EXPECT_DOUBLE_EQ(d(ac::kDa), -m.data.aileron.rate_limit);
  EXPECT_DOUBLE_EQ(d(ac::kDr), m.data.rudder.rate_limit);
  for (int i = 0; i < 400; ++i) {
    const auto dx = m.deriv(0, x, u);
    x.tail<3>() += 0.01 * dx.tail<3>();
  }
  // 4 s exceeds 60 actuator time constants; 1e-12 rad covers accumulated roundoff.
  EXPECT_NEAR(x(ac::kDe), m.data.elevator.pos_max, 1e-12);
  EXPECT_NEAR(x(ac::kDa), m.data.aileron.pos_min, 1e-12);
  EXPECT_NEAR(x(ac::kDr), m.data.rudder.pos_max, 1e-12);
  EXPECT_EQ(m.evaluate(0, x, u).surfaces.tail<3>(), x.tail<3>());
  m.actuators = false;
  EXPECT_EQ(m.deriv(0, x, u).tail<3>(), Vec3::Zero());
  EXPECT_EQ(m.evaluate(0, x, u).surfaces.tail<3>(), Vec3(u.de, u.da, u.dr));
}

TEST(F16Model, ZeroAirspeedDoesNotNaN) {
  auto m = model();
  m.wind.constant_N = Vec3(10, -5, 0);
  auto s = flight_state();
  s(0) = 0;
  const auto x = x17_from_x13(m, s, {});
  const auto e = m.evaluate(0, x, {});
  EXPECT_EQ(e.V, 0);
  EXPECT_EQ(e.alpha, 0);
  EXPECT_EQ(e.beta, 0);
  EXPECT_EQ(e.qbar, 0);
  EXPECT_TRUE(m.deriv(0, x, {}).allFinite());
  EXPECT_TRUE(f13(m, s, {}).allFinite());
}

TEST(F16Model, MonitorFlagsOutOfTableAlpha) {
  auto m = model();
  LookupMonitor monitor;
  m.monitor = &monitor;
  auto s = flight_state();
  s(1) = 50 * units::deg;
  const auto x = x17_from_x13(m, s, {});
  EXPECT_TRUE(m.deriv(0, x, {}).allFinite());
  EXPECT_TRUE(monitor.out_of_domain);
  EXPECT_FALSE(monitor.first_violation.empty());
}

TEST(F16Model, CgTransferUsesTotalDampedForces) {
  auto m = model();
  auto s = flight_state();
  s(2) = -0.1;
  s.segment<3>(6) = Vec3(0.2, 0.3, -0.4);
  const auto x = x17_from_x13(m, s, {});
  const auto before = m.evaluate(0, x, {});
  m.xcg += 0.05;
  const auto after = m.evaluate(0, x, {});
  // CG transfer of total force: pitch=-0.05*cbar*Fz, yaw=+0.05*cbar*Fy; 1e-10 Nm roundoff.
  EXPECT_NEAR(after.fm.M_B.y() - before.fm.M_B.y(), -0.05 * m.data.cbar * before.fm.F_B.z(), 1e-10);
  EXPECT_NEAR(after.fm.M_B.z() - before.fm.M_B.z(), 0.05 * m.data.cbar * before.fm.F_B.y(), 1e-10);
}

TEST(F16Model, MassInertiaAndRotorConventions) {
  auto m = model();
  const auto mp = m.mass_props();
  EXPECT_EQ(mp.mass, m.data.mass);
  EXPECT_EQ(mp.inertia(0, 2), -m.data.Ixz);
  EXPECT_EQ(mp.inertia(2, 0), -m.data.Ixz);
  EXPECT_EQ(mp.inertia.diagonal(), Vec3(m.data.Ixx, m.data.Iyy, m.data.Izz));
  EXPECT_EQ(mp.h_rotor, Vec3(m.data.hx, 0, 0));
  m.include_rotor = false;
  EXPECT_EQ(m.mass_props().h_rotor, Vec3::Zero());
}
