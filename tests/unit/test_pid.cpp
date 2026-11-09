#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include "flightsim/control/pid.hpp"

using flightsim::Pid;
using flightsim::PidGains;

TEST(Pid, ProportionalOnly) {
  Pid p({2, 0, 0, .1, .2, -100, 100});
  EXPECT_DOUBLE_EQ(p.update(3, 1, .02), 4);
}
TEST(Pid, IntegralRamp) {
  Pid p({0, 2, 0, .1, .2, -100, 100});
  // 1e-14 covers roundoff in ten forward-Euler additions of 0.2.
  for (int i = 1; i <= 10; ++i) {
    EXPECT_NEAR(p.update(1, 0, .1), .2 * i, 1e-14);
  }
}
TEST(Pid, DerivativeOnMeasurementNoKick) {
  Pid p({0, 0, 2, .1, .2, -100, 100});
  EXPECT_DOUBLE_EQ(p.update(0, 0, .1), 0);
  EXPECT_DOUBLE_EQ(p.update(10, 0, .1), 0);
  EXPECT_DOUBLE_EQ(p.update(10, 1, .1), -10);
  EXPECT_DOUBLE_EQ(p.update(10, 1, .1), -5);
}
TEST(Pid, SaturationClamps) {
  Pid p({2, 0, 0, .1, .2, -1, 1});
  EXPECT_DOUBLE_EQ(p.update(10, 0, .02), 1);
  EXPECT_DOUBLE_EQ(p.update(-10, 0, .02), -1);
}
TEST(Pid, BackCalculationLimitsWindup) {
  Pid tracked({0, 1, 0, 0, .1, -1, 1});
  Pid untracked({0, 1, 0, 0, .1, -1, 1});
  untracked.antiwindup = false;
  for (int i = 0; i < 5000; ++i) {
    tracked.update(10, 0, .02);
    untracked.update(10, 0, .02);
  }
  double recovery = 0;
  for (int i = 0; i < 10; ++i) {
    recovery = tracked.update(-10, 0, .02);
    EXPECT_DOUBLE_EQ(untracked.update(-10, 0, .02), 1);
  }
  EXPECT_LT(recovery, 0);
}
TEST(Pid, BumplessReset) {
  Pid p({2, 1, 3, .1, .2, -10, 10});
  p.update(3, 1, .02);
  p.reset(.4);
  EXPECT_DOUBLE_EQ(p.update(8, 2, .02), .4);
  p.reset(-.3);
  EXPECT_DOUBLE_EQ(p.update(-2, 7, .02), -.3);
  p.reset(20);
  EXPECT_DOUBLE_EQ(p.update(0, 7, .02), 10);
}
TEST(Pid, RejectsInvalidInputs) {
  EXPECT_THROW((Pid({1, 0, 0, -1, .1, -1, 1})), std::invalid_argument);
  EXPECT_THROW((Pid({1, 0, 0, 0, 0, -1, 1})), std::invalid_argument);
  EXPECT_THROW((Pid({1, 0, 0, 0, .1, 1, -1})), std::invalid_argument);
  Pid p({1, 0, 0, 0, .1, -1, 1});
  EXPECT_THROW(p.update(0, 0, 0), std::invalid_argument);
  EXPECT_THROW(p.update(0, std::numeric_limits<double>::quiet_NaN(), .02), std::invalid_argument);
  EXPECT_THROW(p.reset(std::numeric_limits<double>::infinity()), std::invalid_argument);
}
