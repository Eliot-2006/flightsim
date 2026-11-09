#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "flightsim/control/metrics.hpp"

using namespace flightsim;
namespace {
struct Series {
  std::vector<double> t, y;
};
// Samples y(t) on a 0.01 s grid from 0 to t_end.
template <class F>
Series sample(double t_end, F&& f) {
  Series s;
  for (int k = 0; k <= std::lround(t_end / 0.01); ++k) {
    s.t.push_back(0.01 * k);
    s.y.push_back(f(0.01 * k));
  }
  return s;
}
constexpr double kDt = 0.01;  // grid resolution: times are reported at sample instants
}  // namespace

TEST(Metrics, FirstOrderStepHasNoOvershootAndKnownSettlingTime) {
  // Command at t=2 s, 0 -> 1, time constant 5 s: the 2 % band is entered at 5 ln 50 s after the command.
  const auto s = sample(60, [](double t) { return t < 2 ? 0.0 : 1 - std::exp(-(t - 2) / 5); });
  const auto m = step_metrics(s.t, s.y, 2.0, 0.0, 1.0);
  EXPECT_DOUBLE_EQ(m.step, 1.0);
  EXPECT_EQ(m.overshoot_pct, 0.0);
  EXPECT_NEAR(m.settling_2pct_s, 5 * std::log(50.0), kDt);
  EXPECT_NEAR(m.steady_state_error, std::exp(-(60.0 - 10.0 - 2) / 5), 1e-12);  // max over the last 10 s
}

TEST(Metrics, UnderdampedSecondOrderOvershootMatchesClosedForm) {
  const double zeta = 0.5, wd = std::sqrt(1 - zeta * zeta);
  const auto response = [&](double t) {
    return 1 - std::exp(-zeta * t) * (std::cos(wd * t) + zeta / wd * std::sin(wd * t));
  };
  const auto up = sample(40, response);
  const auto m = step_metrics(up.t, up.y, 0.0, 0.0, 1.0);
  const double expected = 100 * std::exp(-std::numbers::pi * zeta / wd);  // 16.30 %
  // The sampled peak sits within kDt of the analytic one where y is flat to second order: error ~ (wd dt)^2.
  EXPECT_NEAR(m.overshoot_pct, expected, 1e-3);
  // Mirror image (negative step, 1000 -> 900): same overshoot, sign handled.
  const auto down = sample(40, [&](double t) { return 1000 - 100 * response(t); });
  const auto md = step_metrics(down.t, down.y, 0.0, 1000.0, 900.0);
  EXPECT_DOUBLE_EQ(md.step, -100.0);
  EXPECT_NEAR(md.overshoot_pct, expected, 1e-3);
  EXPECT_NEAR(md.settling_2pct_s, m.settling_2pct_s, kDt);
}

TEST(Metrics, NeverSettledIsNaNAndUndershootIsNotOvershoot) {
  const auto s =
      sample(10, [](double t) { return 0.9 * t / 10; });      // still ramping at the end, never above 1
  const auto m = step_metrics(s.t, s.y, 0.0, 0.0, 1.0, 0.0);  // zero-length window: the last sample only
  EXPECT_TRUE(std::isnan(m.settling_2pct_s));
  EXPECT_EQ(m.overshoot_pct, 0.0);
  EXPECT_NEAR(m.steady_state_error, 0.1, 1e-12);
}

TEST(Metrics, StepMetricsRejectZeroStep) {
  const auto s = sample(1, [](double) { return 0.0; });
  EXPECT_THROW(step_metrics(s.t, s.y, 0.0, 1.0, 1.0), std::invalid_argument);
}

TEST(Metrics, DisturbanceDecayRecoversToTwoPercentOfPeak) {
  // Deviation A e^{-(t-3)/4} from t=3 s: peak A at t=3, back inside 2 % of the peak after 4 ln 50 s.
  const double A = 7.0;
  const auto s = sample(80, [&](double t) { return 100 + (t < 3 ? 0.0 : A * std::exp(-(t - 3) / 4)); });
  const auto d = disturbance_metrics(s.t, s.y, 100.0, 3.0, 80.0);
  EXPECT_NEAR(d.peak_deviation, A, 1e-12);
  EXPECT_NEAR(d.recovery_2pct_s, 4 * std::log(50.0), kDt);
  // A window that ends before recovery reports NaN, not the window length.
  EXPECT_TRUE(std::isnan(disturbance_metrics(s.t, s.y, 100.0, 3.0, 10.0).recovery_2pct_s));
  // A zero disturbance recovers instantly.
  const auto flat = sample(5, [](double) { return 100.0; });
  EXPECT_EQ(disturbance_metrics(flat.t, flat.y, 100.0, 0.0, 5.0).recovery_2pct_s, 0.0);
}

TEST(Metrics, PeakAbsDeviationAndRateRespectTheWindow) {
  const auto s = sample(10, [](double t) { return t < 5 ? 3 * t : 15 - 8 * (t - 5); });  // peak 15 at t=5
  EXPECT_NEAR(peak_abs_deviation(s.t, s.y, 0.0, 0.0, 10.0), 25.0, 1e-9);                 // at t=10: 15-40
  EXPECT_NEAR(peak_abs_deviation(s.t, s.y, 0.0, 0.0, 5.0), 15.0, 1e-9);
  EXPECT_NEAR(peak_abs_rate(s.t, s.y), 8.0, 1e-9);
}
