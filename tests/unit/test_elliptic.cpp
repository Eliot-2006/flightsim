#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "flightsim/math/elliptic.hpp"

using namespace flightsim;

TEST(Elliptic, MatchesScipy) {
  std::ifstream stream(std::string(FLIGHTSIM_SOURCE_DIR) + "/tests/data/ellipj_reference.csv");
  ASSERT_TRUE(stream);
  std::string line;
  std::getline(stream, line);
  int count = 0;
  while (std::getline(stream, line)) {
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream row(line);
    double u, m, sn, cn, dn;
    ASSERT_TRUE(static_cast<bool>(row >> u >> m >> sn >> cn >> dn));
    const auto actual = jacobi_elliptic(u, m);
    // Task 2.2: independent SciPy double-precision reference, 1e-12 absolute bound.
    EXPECT_NEAR(actual.sn, sn, 1e-12) << u << ", " << m;
    EXPECT_NEAR(actual.cn, cn, 1e-12) << u << ", " << m;
    EXPECT_NEAR(actual.dn, dn, 1e-12) << u << ", " << m;
    ++count;
  }
  EXPECT_EQ(count, 606);
}

TEST(Elliptic, Identities) {
  for (double m : {0., .1, .5, .9, .99, .999, 1. - 1e-12}) {
    for (int i = -100; i <= 100; ++i) {
      const auto [sn, cn, dn] = jacobi_elliptic(i * .1, m);
      // Task 2.2: Jacobi algebraic identities to 1e-14 (~45 double ulp).
      EXPECT_NEAR(sn * sn + cn * cn, 1, 1e-14);
      EXPECT_NEAR(dn * dn + m * sn * sn, 1, 1e-14);
    }
  }
}

TEST(Elliptic, LimitsM0AndM1) {
  for (double u : {-5., -1., 0., 1., 5.}) {
    const auto zero = jacobi_elliptic(u, 0);
    // Task 2.2 m=0 trigonometric limit, within a few ulp.
    EXPECT_NEAR(zero.sn, std::sin(u), 1e-15);
    EXPECT_NEAR(zero.cn, std::cos(u), 1e-15);
    EXPECT_EQ(zero.dn, 1);
    const auto one = jacobi_elliptic(u, 1. - 1e-12);
    // m->1 differs from tanh/sech by O((1-m)*exp(|u|)); bound 1e-10 on |u|<=5.
    EXPECT_NEAR(one.sn, std::tanh(u), 1e-10);
    EXPECT_NEAR(one.cn, 1 / std::cosh(u), 1e-10);
    EXPECT_NEAR(one.dn, 1 / std::cosh(u), 1e-10);
  }
}

TEST(Elliptic, RejectsInvalidArguments) {
  for (double m : {-0.1, 1., 2., std::nan("")}) EXPECT_THROW(jacobi_elliptic(1, m), std::invalid_argument);
  EXPECT_THROW(jacobi_elliptic(std::numeric_limits<double>::infinity(), .5), std::invalid_argument);
}
