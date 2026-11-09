#include <gtest/gtest.h>

#include <cmath>

#include "flightsim/math/frames.hpp"
#include "flightsim/math/quat.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;

TEST(Frames, AlphaPositiveForPositiveW) {
  const AirData a = air_data_from_body_velocity(Vec3(100.0, 0.0, 10.0));
  EXPECT_GT(a.alpha, 0.0);
  EXPECT_NEAR(a.alpha, std::atan2(10.0, 100.0), 1e-15);  // definition alpha = atan2(w, u) (PRD 4.1)
  EXPECT_EQ(a.beta, 0.0);
  EXPECT_NEAR(a.V, std::sqrt(100.0 * 100.0 + 10.0 * 10.0), 1e-13);
}

TEST(Frames, BetaPositiveForPositiveV) {
  const AirData a = air_data_from_body_velocity(Vec3(100.0, 10.0, 0.0));
  EXPECT_GT(a.beta, 0.0);
  EXPECT_NEAR(a.beta, std::asin(10.0 / std::sqrt(10100.0)), 1e-15);  // definition beta = asin(v/V) (PRD 4.1)
  EXPECT_EQ(a.alpha, 0.0);
}

TEST(Frames, ZeroAirspeedDefined) {
  for (const Vec3& v : {Vec3(0, 0, 0), Vec3(1e-12, -1e-12, 1e-12)}) {
    const AirData a = air_data_from_body_velocity(v);
    EXPECT_TRUE(std::isfinite(a.V) && std::isfinite(a.alpha) && std::isfinite(a.beta));
    EXPECT_EQ(a.alpha, 0.0);  // V < 1e-9 -> alpha = beta = 0 (documented in the interface)
    EXPECT_EQ(a.beta, 0.0);
  }
  EXPECT_EQ(air_data_from_body_velocity(Vec3::Zero()).V, 0.0);
}

TEST(Frames, EulerRateMatrixMatchesFiniteDifference) {
  const double phi = 0.3, theta = 0.5, psi = -0.7;
  const Vec3 w(0.2, -0.4, 0.3);
  const Vec4 q = quat_from_euler(phi, theta, psi);
  const Vec4 qdot = 0.5 * quat_mul(q, Vec4(0.0, w.x(), w.y(), w.z()));  // PRD 4.2 kinematics
  const double d = 1e-6;
  const Vec3 fd =
      (euler_from_quat(quat_normalized(q + d * qdot)) - euler_from_quat(quat_normalized(q - d * qdot))) /
      (2 * d);
  // Tolerance 1e-7 (plan value): central-difference truncation ~d^2 = 1e-12 plus roundoff ~1e-16/d = 1e-10.
  EXPECT_LE((euler_rate_matrix(phi, theta) * w - fd).norm(), 1e-7);
}

TEST(Units, ConstantsAreConsistent) {
  EXPECT_EQ(units::ft, 0.3048);   // exact by definition
  EXPECT_EQ(units::g0, 9.80665);  // exact by definition (standard gravity)
  EXPECT_EQ(units::deg * 180.0, std::numbers::pi);
  // slug = lbf / (ft/s^2). The PRD gives slug to 10 significant digits, so the relative mismatch is below
  // 5e-10; 1e-9 is that rounding with ~2x margin.
  EXPECT_NEAR(units::slug / (units::lbf / units::ft), 1.0, 1e-9);
}
