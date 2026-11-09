#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>

#include "flightsim/analysis/flight_vars.hpp"

using namespace flightsim;

namespace {
F16Model model() {
  F16Model m{load_f16_data(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml")};
  m.grav.g0 = m.data.g_sl;
  m.wind.constant_N = Vec3(10, -5, 2);
  return m;
}

Vec<13> flight_state() {
  Vec<13> x;
  x << 152.4, 0.13, -0.08, -0.4, 0.3, 0.7, 0.2, -0.15, 0.1, 100, -200, 3048, 60;
  return x;
}
}  // namespace

TEST(FlightVars, RoundTripX13) {
  auto m = model();
  const F16Input u{0.7, -0.03, 0.02, -0.01};
  const auto x13 = flight_state();
  for (bool actuators : {false, true}) {
    m.actuators = actuators;
    const auto x17 = x17_from_x13(m, x13, u);
    const auto roundtrip = x13_from_x17(m, x17);
    // 1e-12 relative/scale-one floor: roundoff in trig/DCM/inverse trig (Task 4.5).
    for (int i = 0; i < 13; ++i) {
      EXPECT_NEAR(roundtrip(i), x13(i), 1e-12 * std::max(1.0, std::abs(x13(i)))) << i;
    }
    EXPECT_DOUBLE_EQ(x17(ac::kDe), u.de);
    EXPECT_DOUBLE_EQ(x17(ac::kDa), u.da);
    EXPECT_DOUBLE_EQ(x17(ac::kDr), u.dr);
  }
}

TEST(FlightVars, F13EulerRatesMatchFiniteDifference) {
  const auto m = model();
  const F16Input u{0.7, -0.03, 0.02, -0.01};
  const auto x13 = flight_state();
  const auto x17 = x17_from_x13(m, x13, u);
  const auto d = m.deriv(0, x17, u);
  // Central FD at epsilon=1e-5 s: O(epsilon^2) truncation plus roundoff in altitude subtraction.
  constexpr double eps = 1e-5;
  const Vec<13> fd =
      (x13_from_x17(m, Vec<17>(x17 + eps * d)) - x13_from_x17(m, Vec<17>(x17 - eps * d))) / (2 * eps);
  const auto actual = f13(m, x13, u);
  for (int i = 0; i < 13; ++i) {
    // Task 4.5: 1e-7 absolute; tests all components, including rotating-frame air acceleration.
    EXPECT_NEAR(actual(i), fd(i), 1e-7) << i;
  }
}
