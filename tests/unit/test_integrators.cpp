#include <gtest/gtest.h>

#include <cmath>

#include "flightsim/core/integrators.hpp"
#include "flightsim/core/state.hpp"

using namespace flightsim;

namespace {
Vec<1> decay(double, const Vec<1>& x) { return -x; }
}  // namespace

TEST(Integrators, EulerScalarExponential) {
  const Vec<1> x1 = euler_step(decay, 0.0, Vec<1>(1.0), 0.1);
  EXPECT_NEAR(x1(0), 0.9, 1e-16);  // 1 + h*(-1) = 0.9 up to the rounding of 0.1
}

TEST(Integrators, Rk4ScalarExponential) {
  const double h = 0.1;
  const Vec<1> x1 = rk4_step(decay, 0.0, Vec<1>(1.0), h);
  // One-step truncation error ~ h^5/120 = 8.3e-8 -> bound 1e-7 (plan value).
  EXPECT_LE(std::abs(x1(0) - std::exp(-h)), 1e-7);
  // RK4 on x' = -x is exactly the degree-4 Taylor polynomial; 1e-15 = a few ulp of 0.9.
  EXPECT_NEAR(x1(0), 1 - h + h * h / 2 - h * h * h / 6 + h * h * h * h / 24, 1e-15);
}

TEST(Integrators, RenormalizeLeavesOthers) {
  Vec<14> x = Vec<14>::Constant(7.0);
  rb::quat(x) = Vec4(2, 0, 0, 0);
  const Vec<14> before = x;
  renormalize_quat(x);
  EXPECT_EQ(rb::quat(x), Vec4(1, 0, 0, 0));  // exact: 2/2
  for (int i = 0; i < 14; ++i)
    if (i < rb::kQuat || i >= rb::kQuat + 4) EXPECT_EQ(x(i), before(i)) << i;
}

TEST(Integrators, StepDispatcherCountsEvaluations) {
  StepStats st;
  Vec<13> x = Vec<13>::Zero();
  rb::quat(x) = Vec4(1, 0, 0, 0);
  auto f = [](double, const Vec<13>&) { return Vec<13>(Vec<13>::Zero()); };
  step(Integrator::Euler, f, 0.0, x, 0.1, true, st);
  EXPECT_EQ(st.nfev, 1);
  step(Integrator::RK4, f, 0.0, x, 0.1, true, st);
  EXPECT_EQ(st.nfev, 5);
  step(Integrator::RKMK4, f, 0.0, x, 0.1, true, st);
  EXPECT_EQ(st.nfev, 9);
  step(Integrator::LieEuler, f, 0.0, x, 0.1, true, st);
  EXPECT_EQ(st.nfev, 10);
  step(Integrator::SymplecticEuler, f, 0.0, x, 0.1, true, st);
  EXPECT_EQ(st.nfev, 11);
}

TEST(Integrators, NamesRoundTrip) {
  for (Integrator i : {Integrator::Euler, Integrator::RK4, Integrator::SymplecticEuler, Integrator::LieEuler,
                       Integrator::RKMK4})
    EXPECT_EQ(integrator_from_string(to_string(i)), i);
  EXPECT_THROW(integrator_from_string("bogus"), std::invalid_argument);
}
