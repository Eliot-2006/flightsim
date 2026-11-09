#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include "flightsim/env/wind.hpp"

using namespace flightsim;

TEST(Wind, ConstantOnly) {
  WindModel wind{Vec3(10, -5, 2), std::nullopt};
  EXPECT_EQ(wind.wind_N(-10), wind.constant_N);
  EXPECT_EQ(wind.wind_N(100), wind.constant_N);
}

TEST(Wind, GustZeroBeforeAndAfter) {
  WindModel wind{Vec3(1, 2, 3), Gust{2, 4, 10, Vec3::UnitZ()}};
  EXPECT_EQ(wind.wind_N(1), wind.constant_N);
  EXPECT_EQ(wind.wind_N(7), wind.constant_N);
}

TEST(Wind, GustPeakAtMidpoint) {
  WindModel wind{Vec3::Zero(), Gust{2, 4, 10, Vec3::UnitY()}};
  // Task 4.4: peak amplitude within 1e-15 m/s (cos(pi) is exactly -1 in double).
  EXPECT_LE((wind.wind_N(4) - Vec3(0, 10, 0)).norm(), 1e-15);
}

TEST(Wind, GustContinuous) {
  WindModel wind{Vec3::Zero(), Gust{2, 4, 10, Vec3::UnitZ()}};
  EXPECT_EQ(wind.wind_N(2), Vec3::Zero());
  EXPECT_EQ(wind.wind_N(6), Vec3::Zero());
  // Taylor bound A*pi^2*(epsilon/T)^2 for epsilon=1e-6 s, with rounding slack.
  EXPECT_LT(wind.wind_N(2 + 1e-6).norm(), 7e-12);
  EXPECT_LT(wind.wind_N(6 - 1e-6).norm(), 7e-12);
}

TEST(Wind, DirectionNormalized) {
  WindModel wind{Vec3::Zero(), Gust{0, 2, 5, Vec3(0, 3, 4)}};
  // Exact 3-4-5 direction; 1e-15 m/s covers normalization arithmetic.
  EXPECT_LE((wind.wind_N(1) - Vec3(0, 3, 4)).norm(), 1e-15);
}

TEST(Wind, RejectsInvalidGust) {
  WindModel wind{Vec3::Zero(), Gust{0, 0, 5, Vec3::UnitX()}};
  EXPECT_THROW(wind.wind_N(0), std::invalid_argument);
  wind.gust->duration_s = -1;
  EXPECT_THROW(wind.wind_N(0), std::invalid_argument);
  wind.gust->duration_s = 1;
  wind.gust->direction_N.setZero();
  EXPECT_THROW(wind.wind_N(0), std::invalid_argument);
  wind.gust->direction_N = Vec3::UnitX();
  wind.gust->amplitude_mps = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(wind.wind_N(0), std::invalid_argument);
}
