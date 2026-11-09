#include <gtest/gtest.h>

#include "flightsim/analysis/jacobian.hpp"

using namespace flightsim;

TEST(Jacobian, JacobianOfLinearMapExact) {
  Eigen::Matrix<double, 2, 3> A;
  A << 1, 2, -3, 4, -5, 6;
  auto f = [&](const Vec3& x) -> Vec<2> { return A * x + Vec<2>(2, -1); };
  const auto J = jacobian_central<2, 3>(f, Vec3(0, 2, -3), Vec3(1, 1, 1));
  // Task 2.3: central differences of a linear map have only subtraction roundoff, <=1e-9.
  EXPECT_LE((J - A).norm(), 1e-9);
}

TEST(Jacobian, ScalesZeroComponentsAndRejectsInvalidSteps) {
  auto f = [](const Vec3& x) -> Vec3 { return x.array().sin().matrix(); };
  // At x=0, scales supply the nonzero stencil; sin'(0)=1, O(eps^2) truncation <1e-9.
  EXPECT_LE((jacobian_central<3, 3>(f, Vec3::Zero(), Vec3(1, 2, 3)) - Mat3::Identity()).norm(), 1e-9);
  const Vec3 x = Vec3::Zero(), scale = Vec3::Ones();
  EXPECT_THROW((jacobian_central<3, 3>(f, x, scale, 0)), std::invalid_argument);
  EXPECT_THROW((jacobian_central<3, 3>(f, x, Vec3::Zero())), std::invalid_argument);
}
