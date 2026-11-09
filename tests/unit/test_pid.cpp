#include <gtest/gtest.h>

#include <algorithm>
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
TEST(Pid, PerSampleLimitsAreTrackedLikeConstructedLimits) {
  // A stage downstream of the PID (the autopilot's combined theta bound) can be tighter than the gains'
  // limits; back-calculation must then track that tighter bound exactly as if it were the gains' own.
  Pid tight({0.5, 1, 0, 0, .1, -1, 1});
  Pid wide({0.5, 1, 0, 0, .1, -10, 10});
  for (int i = 0; i < 500; ++i) {
    EXPECT_DOUBLE_EQ(wide.update(10, 0, .02, -1, 1), tight.update(10, 0, .02)) << i;
  }
  for (int i = 0; i < 50; ++i) {
    EXPECT_DOUBLE_EQ(wide.update(-10, 0, .02, -1, 1), tight.update(-10, 0, .02)) << i;
  }
  EXPECT_THROW(wide.update(0, 0, .02, 1, -1), std::invalid_argument);
}
TEST(Pid, ProportionalOnlyLoopKeepsNoBiasAfterSaturationOrReset) {
  // ki == 0 (the autopilot's pitch loop) has no integrator: neither back-calculation while saturated nor a
  // bumpless reset may store an offset, which nothing would ever integrate away again.
  Pid p({2, 0, 0, 0, .5, -1, 1});
  for (int i = 0; i < 500; ++i) {
    EXPECT_DOUBLE_EQ(p.update(10, 0, .02), 1);  // 10 s saturated
  }
  EXPECT_DOUBLE_EQ(p.update(.2, 0, .02), .4);  // out of saturation: exactly kp e
  p.reset(.7);
  EXPECT_DOUBLE_EQ(p.update(.1, 0, .02), .7);  // the reset output, for that sample only
  EXPECT_DOUBLE_EQ(p.update(.1, 0, .02), .2);  // then kp e again, no stored bias
}
TEST(Pid, BackCalculationOnAppliedOutputBoundsIntegratorUnderRateLimit) {
  // A rate-limited actuator (0.5/s) lags a PI whose command never reaches its own +-10 limits. Tracking the
  // command saturation sees nothing to correct, so that integrator winds up as ki e t. Tracking the applied
  // (actuator) output holds the command within ki e tt (= 0.1) of the actuator.
  const PidGains g{0, 1, 0, 0, .1, -10, 10};
  Pid command(g), applied(g);
  double actuator = 0, c = 0, u = 0;
  for (int i = 0; i < 250; ++i) {  // 5 s at 50 Hz, constant unit error
    c = command.update(1, 0, .02);
    u = applied.update(1, 0, .02, -10, 10, actuator);
    actuator += std::clamp(u - actuator, -.5 * .02, .5 * .02);
  }
  EXPECT_NEAR(actuator, 2.5, 1e-12);  // rate-limited the whole time
  EXPECT_NEAR(c, 5, 1e-12);           // wound up: ki e t
  EXPECT_GT(u - actuator, 0);
  EXPECT_LT(u - actuator, .1);  // bounded by ki e tt
  // Reversing the error takes the applied-tracking command below the actuator within 0.2 s; the wound-up
  // one is still 2 above it.
  for (int i = 0; i < 10; ++i) {
    c = command.update(-1, 0, .02);
    u = applied.update(-1, 0, .02, -10, 10, actuator);
    actuator += std::clamp(u - actuator, -.5 * .02, .5 * .02);
  }
  EXPECT_LT(u, actuator);
  EXPECT_GT(c, actuator + 2);
  EXPECT_THROW(applied.update(0, 0, .02, -1, 1, std::numeric_limits<double>::quiet_NaN()),
               std::invalid_argument);
}
