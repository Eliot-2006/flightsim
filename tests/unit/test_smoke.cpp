#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "flightsim/version.hpp"

TEST(Smoke, VersionStringNonEmpty) { EXPECT_FALSE(flightsim::version_string().empty()); }

TEST(Smoke, EigenUnitXDotSelfIsOne) {
  // Exact: 1*1 + 0*0 + 0*0 has no rounding.
  EXPECT_EQ(Eigen::Vector3d::UnitX().dot(Eigen::Vector3d::UnitX()), 1.0);
}
