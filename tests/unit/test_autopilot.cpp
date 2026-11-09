#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "flightsim/control/autopilot.hpp"

using namespace flightsim;
namespace {
const std::filesystem::path root = FLIGHTSIM_SOURCE_DIR;
AutopilotGains committed() { return load_autopilot_gains(root / "scenarios/control/autopilot_gains.toml"); }

// Trim at the design point (values from docs/design.md "M7a: control operating point"): theta = alpha =
// 3.367607324 deg, elevator -0.656145669 deg, throttle 0.168247601, V = 153.0096 m/s, h = 3048 m.
constexpr double kThetaTrim = 3.367607324 * 3.14159265358979323846 / 180;
constexpr double kDeTrim = -0.656145669 * 3.14159265358979323846 / 180;
constexpr double kThrottleTrim = 0.168247601;
const F16Input kTrim{kThrottleTrim, kDeTrim, 0, 0};
const ApMeasurement kSteady{3048, 153.0096, kThetaTrim, 0, kDeTrim};  // elevator resting at trim
const ApCommand kHold{3048, 153.0096};
}  // namespace

TEST(Autopilot, LoadsCommittedGains) {
  const auto g = committed();
  EXPECT_DOUBLE_EQ(g.rate_hz, 50.0);
  EXPECT_DOUBLE_EQ(g.pitch_rate.kp, -1.174838697);   // negative: S&L positive elevator is nose down
  EXPECT_DOUBLE_EQ(g.pitch_rate.umax, 0.436332313);  // +-25 deg elevator position limit
  EXPECT_DOUBLE_EQ(g.pitch.ki, 0.0);
  EXPECT_DOUBLE_EQ(g.altitude.umax, 0.2617993878);  // +-15 deg theta_cmd
  EXPECT_DOUBLE_EQ(g.airspeed.umin, 0.0);
  EXPECT_DOUBLE_EQ(g.airspeed.umax, 1.0);
  EXPECT_DOUBLE_EQ(g.airspeed_to_pitch, -0.001094945016);
  EXPECT_DOUBLE_EQ(g.reference_airspeed_mps, 153.0096);
}

TEST(Autopilot, GainsLoaderRejectsUnknownAndMissingKeys) {
  const auto path = std::filesystem::temp_directory_path() / "flightsim_test_gains.toml";
  std::ifstream in(root / "scenarios/control/autopilot_gains.toml");
  const std::string good((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  const auto write = [&](const std::string& text) { std::ofstream(path) << text; };
  write(good + "\n[extra]\nx = 1\n");
  EXPECT_THROW(load_autopilot_gains(path), std::runtime_error);
  write(good.substr(0, good.find("[crossfeed]")));  // crossfeed table missing
  EXPECT_THROW(load_autopilot_gains(path), std::runtime_error);
  const auto kd = good.find("kd = 0.0");
  write(good.substr(0, kd) + "kdd = 0.0" + good.substr(kd + 8));  // misspelled key
  EXPECT_THROW(load_autopilot_gains(path), std::runtime_error);
}

TEST(Autopilot, EngageIsBumplessAndHoldsTrimAtZeroError) {
  Autopilot ap(committed(), true);
  ap.engage(kSteady, kTrim);
  for (int i = 0; i < 500; ++i) {
    const auto o = ap.update(kSteady, kHold);
    // Absolute commands: reset to the trim outputs, with no second trim offset added on top.
    EXPECT_DOUBLE_EQ(o.u.throttle, kThrottleTrim) << i;
    EXPECT_DOUBLE_EQ(o.u.de, kDeTrim) << i;
    EXPECT_DOUBLE_EQ(o.theta_cmd, kThetaTrim) << i;
    EXPECT_DOUBLE_EQ(o.q_cmd, 0.0) << i;
    EXPECT_EQ(o.u.da, 0.0);
    EXPECT_EQ(o.u.dr, 0.0);
  }
}

TEST(Autopilot, CrossfeedIsMeasuredAirspeedTimesSlopeOnThePitchReference) {
  const auto g = committed();
  Autopilot ap(g, true);
  ap.engage(kSteady, kTrim);
  ap.update(kSteady, kHold);
  ApMeasurement fast = kSteady;
  fast.V += 10;  // measured V, not the commanded one: the command stays at V_trim
  const auto o = ap.update(fast, kHold);
  // Altitude error is zero, so theta_alt stays at trim theta and theta_cmd moves by exactly k_V * dV.
  EXPECT_NEAR(o.theta_cmd - kThetaTrim, g.airspeed_to_pitch * 10, 1e-15);
  EXPECT_LT(o.theta_cmd, kThetaTrim);  // faster -> lower trim pitch (design.md: k_V < 0)
}

TEST(Autopilot, CascadeSignsAndGainsForAnAltitudeStepAreHandComputable) {
  const auto g = committed();
  const double T = 1 / g.rate_hz;
  Autopilot ap(g, true);
  ap.engage(kSteady, kTrim);
  ap.update(kSteady, kHold);  // sample 0: the reset outputs
  const double dh = 30.48;
  const auto o = ap.update(kSteady, {3048 + dh, 153.0096});
  // PI with forward-Euler integral applied in the same sample (docs/design.md "M7a: PID discrete form").
  const double theta_cmd = kThetaTrim + (g.altitude.kp + T * g.altitude.ki) * dh;
  EXPECT_NEAR(o.theta_cmd, theta_cmd, 1e-15);
  EXPECT_NEAR(o.q_cmd, g.pitch.kp * (theta_cmd - kThetaTrim), 1e-15);
  EXPECT_NEAR(o.u.de, kDeTrim + (g.pitch_rate.kp + T * g.pitch_rate.ki) * o.q_cmd, 1e-15);
  EXPECT_GT(o.q_cmd, 0);                        // climb -> pitch up
  EXPECT_LT(o.u.de, kDeTrim);                   // trailing-edge-up elevator (S&L: positive = down)
  EXPECT_NEAR(o.u.throttle, kThrottleTrim, 0);  // airspeed loop sees no error
}

TEST(Autopilot, AirspeedStepRaisesThrottleOnly) {
  const auto g = committed();
  Autopilot ap(g, true);
  ap.engage(kSteady, kTrim);
  ap.update(kSteady, kHold);
  const double dv = 5;  // small enough not to saturate the throttle (kp dV = 0.36)
  const auto o = ap.update(kSteady, {3048, 153.0096 + dv});
  EXPECT_NEAR(o.u.throttle, kThrottleTrim + (g.airspeed.kp + g.airspeed.ki / g.rate_hz) * dv, 1e-15);
  EXPECT_DOUBLE_EQ(o.theta_cmd, kThetaTrim);  // altitude loop unaffected; measured V still at V_trim
}

TEST(Autopilot, CombinedThetaBoundIncludesTheCrossfeed) {
  const auto g = committed();
  const double bound = g.altitude.umax;
  for (const double dv : {-40.0, 0.0, 40.0}) {
    for (const double sign : {1.0, -1.0}) {
      Autopilot ap(g, true);
      ap.engage(kSteady, kTrim);
      ApMeasurement m = kSteady;
      m.V += dv;
      for (int i = 0; i < 2000; ++i) {
        const auto o = ap.update(m, {3048 + sign * 1e4, 153.0096});
        // theta_cmd = clamp(theta_altitude_PI + k_V (V - V_trim), -15 deg, +15 deg), also at dV = +-40 m/s
        // where the cross-feed alone is +-2.5 deg and the altitude PI sits on its own bound.
        ASSERT_LE(std::abs(o.theta_cmd), bound) << "dv=" << dv << " i=" << i;
      }
      const auto o = ap.update(m, {3048 + sign * 1e4, 153.0096});
      // The PI sits on its own bound (or on the tighter bound that leaves room for the cross-feed), so the
      // sum is the PI bound plus the cross-feed, limited to +-15 deg.
      const double expected = std::clamp(sign * bound + g.airspeed_to_pitch * dv, -bound, bound);
      EXPECT_NEAR(o.theta_cmd, expected, 1e-15) << "dv=" << dv << " sign=" << sign;
    }
  }
}

TEST(Autopilot, AntiwindupFlagControlsBackCalculation) {
  const auto g = committed();
  const auto saturate_then_release = [&](bool antiwindup) {
    Autopilot ap(g, antiwindup);
    ap.engage(kSteady, kTrim);
    ap.update(kSteady, kHold);        // the first sample after engage outputs the reset values
    for (int i = 0; i < 1000; ++i) {  // 20 s with a 300 m error: theta_cmd pinned at +15 deg
      ap.update(kSteady, {3048 + 300, 153.0096});
    }
    return ap.update(kSteady, kHold).theta_cmd;  // error removed: how much integral state remains
  };
  // Back-calculation drags the integrator down while saturated; without it the integrator keeps growing.
  EXPECT_LT(saturate_then_release(true), saturate_then_release(false));
}

TEST(Autopilot, PitchRateBackCalculationTracksTheMeasuredElevator) {
  // The pitch-rate PI tracks the applied elevator (the actuator state), not its own command. With zero error
  // and the elevator held 5 deg from the command, back-calculation relaxes the command toward it with time
  // constant tt; without anti-windup it stays at trim.
  const auto g = committed();
  const double offset = 5 * 3.14159265358979323846 / 180;
  const auto run = [&](bool antiwindup) {
    Autopilot ap(g, antiwindup);
    ap.engage(kSteady, kTrim);
    ApMeasurement lagging = kSteady;
    lagging.de = kDeTrim + offset;
    double de = 0;
    for (int i = 0; i < 250; ++i) {  // 5 s
      de = ap.update(lagging, kHold).u.de;
    }
    return de;
  };
  // Output n (1-based) is de_trim + offset (1 - (1 - T/tt)^(n-1)): the reset sample outputs trim exactly.
  const double expected = kDeTrim + offset * (1 - std::pow(1 - 1 / (g.rate_hz * g.pitch_rate.tt), 249));
  EXPECT_NEAR(run(true), expected, 1e-12);
  EXPECT_DOUBLE_EQ(run(false), kDeTrim);
}

TEST(Autopilot, UpdateBeforeEngageThrows) {
  Autopilot ap(committed(), true);
  EXPECT_THROW(ap.update(kSteady, kHold), std::logic_error);
}
