#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "flightsim/validation/analytic.hpp"
#include "flightsim/validation/cases.hpp"
#include "flightsim/validation/convergence.hpp"

using namespace flightsim;
using namespace flightsim::validation;

TEST(Rotational, R2ReferenceSatisfiesKinematics) {
  const auto c = case_R2();
  EXPECT_LE((c.exact(0) - c.x0).norm(), 1e-14);  // Closed-form initial condition, ~45 ulp.
  for (double t : {0., .3, 1., 4., 9.}) {
    const Vec<13> x = c.exact(t);
    const Vec4 fd = (rb::quat(c.exact(t + 1e-5)) - rb::quat(c.exact(t - 1e-5))) / 2e-5;
    const Vec<13> d = rigid_body_deriv(x, c.mp, {}, c.g_N);
    // PRD 8.3 reference self-check: central FD truncation O(dt^2), residual <=1e-9.
    EXPECT_LE((fd - rb::quat(d)).norm(), 1e-9) << t;
    const Vec3 wfd = (rb::omega(c.exact(t + 1e-5)) - rb::omega(c.exact(t - 1e-5))) / 2e-5;
    EXPECT_LE((wfd - rb::omega(d)).norm(), 1e-9);
  }
}

TEST(Rotational, R3ReferenceSatisfiesEuler) {
  const auto c = case_R3();
  const Vec3 I(1, 2, 3);
  const Vec3 w0 = rb::omega(c.x0);
  const double E0 = .5 * w0.dot(I.cwiseProduct(w0));
  const double L0 = I.cwiseProduct(w0).norm();
  EXPECT_LE((exact_asymmetric_omega(I, w0, 0) - w0).norm(), 1e-14);  // Initial-condition roundoff.
  EXPECT_GT(exact_asymmetric_omega(I, w0, 1e-3).y(), 0);             // Euler: w2dot(0)=+0.5, not -0.5.
  for (int i = 0; i <= 200; ++i) {
    const double t = i * .1;
    const Vec3 w = exact_asymmetric_omega(I, w0, t);
    const Vec3 fd =
        (exact_asymmetric_omega(I, w0, t + 1e-5) - exact_asymmetric_omega(I, w0, t - 1e-5)) / 2e-5;
    const Vec3 euler = -I.cwiseInverse().cwiseProduct(w.cross(I.cwiseProduct(w)));
    // PRD 8.3: FD residual <=1e-9 and invariant absolute errors <=1e-13.
    EXPECT_LE((fd - euler).norm(), 1e-9) << t;
    EXPECT_NEAR(.5 * w.dot(I.cwiseProduct(w)), E0, 1e-13);
    EXPECT_NEAR(I.cwiseProduct(w).norm(), L0, 1e-13);
  }
  // Initial signs also determine the direction of motion on the same LL branch.
  for (double s1 : {-1., 1.})
    for (double s3 : {-1., 1.}) {
      const Vec3 initial(s1 * .5, 0, s3);
      const Vec3 fd =
          (exact_asymmetric_omega(I, initial, 1e-5) - exact_asymmetric_omega(I, initial, -1e-5)) / 2e-5;
      EXPECT_NEAR(fd.y(), s1 * s3 * .5, 1e-9);  // Same PRD Euler residual bound.
    }
}

namespace {
void expect_order(const std::vector<double>& h, const std::vector<double>& e, double scale,
                  Integrator method) {
  const auto fit = observed_order(h, e, scale);
  // PRD 10.3: first order in [.95,1.05], fourth order in [3.9,4.1].
  const bool fourth = method == Integrator::RK4 || method == Integrator::RKMK4;
  EXPECT_GE(fit.order, fourth ? 3.9 : .95);
  EXPECT_LE(fit.order, fourth ? 4.1 : 1.05);
  std::printf("[order] %s %.8f levels %d-%d\n", to_string(method), fit.order, fit.first_level,
              fit.last_level);
}
void rotational_order(const CaseDef& c, Integrator method) {
  std::vector<double> h, att, omega, momentum;
  const Vec3 H0 = angular_momentum_N(c.x0, c.mp);
  for (int k = 0; k < 8; ++k) {
    // R3 first-order runs start at .025 to exclude the nonlinear pre-asymptotic regime.
    const double h0 =
        c.name == "R3" && (method == Integrator::Euler || method == Integrator::LieEuler) ? .025 : .1;
    h.push_back(h0 / std::ldexp(1., k));
    const auto x = integrate_case(c, method, h.back()).x_final;
    if (c.exact) att.push_back(state_error(x, c.exact(c.t_final)).att);
    const Vec3 wref = c.exact_omega ? c.exact_omega(c.t_final) : Vec3(rb::omega(c.exact(c.t_final)));
    omega.push_back((rb::omega(x) - wref).norm());
    momentum.push_back((angular_momentum_N(x, c.mp) - H0).norm() / H0.norm());
  }
  if (c.name == "R1")
    expect_order(h, att, 1, method);
  else {
    if (c.exact) expect_order(h, att, 1, method);
    expect_order(h, omega, rb::omega(c.x0).norm(), method);
    if (c.name == "R3") expect_order(h, momentum, 1, method);
  }
}
}  // namespace

TEST(Rotational, R1Rk4Order4) { rotational_order(case_R1(), Integrator::RK4); }
TEST(Rotational, R1EulerOrder1) {
  const auto c = case_R1();
  std::vector<double> h, e;
  for (int k = 0; k < 8; ++k) {
    h.push_back(.1 / std::ldexp(1., k));
    // Constant-rate Euler angle superconverges at order 2. The unprojected R^4 norm error is order 1.
    const auto x = integrate_case(c, Integrator::Euler, h.back(), false).x_final;
    e.push_back(std::abs(rb::quat(x).norm() - 1));
  }
  expect_order(h, e, 1, Integrator::Euler);
}
TEST(Rotational, R2Rk4Order4) { rotational_order(case_R2(), Integrator::RK4); }
TEST(Rotational, R2Rkmk4Order4) { rotational_order(case_R2(), Integrator::RKMK4); }
TEST(Rotational, R2LieEulerOrder1) { rotational_order(case_R2(), Integrator::LieEuler); }
TEST(Rotational, R2EulerOrder1) { rotational_order(case_R2(), Integrator::Euler); }
TEST(Rotational, R3Rk4Order4) { rotational_order(case_R3(), Integrator::RK4); }
TEST(Rotational, R3Rkmk4Order4) { rotational_order(case_R3(), Integrator::RKMK4); }
TEST(Rotational, R3LieEulerOrder1) { rotational_order(case_R3(), Integrator::LieEuler); }
TEST(Rotational, R3EulerOrder1) { rotational_order(case_R3(), Integrator::Euler); }

TEST(Rotational, HNConservedRk4) {
  const auto c = case_R3();
  const auto x = integrate_case(c, Integrator::RK4, 1e-3).x_final;
  const Vec3 H0 = angular_momentum_N(c.x0, c.mp);
  // Task 2.2 bound: R3 RK4 at h=.001 conserves inertial momentum to relative 1e-10.
  EXPECT_LE((angular_momentum_N(x, c.mp) - H0).norm() / H0.norm(), 1e-10);
}

TEST(Rotational, R5RotorHNConserved) {
  const auto c = case_R5();
  const auto x = integrate_case(c, Integrator::RK4, 1e-3).x_final;
  const Vec3 H0 = angular_momentum_N(c.x0, c.mp);
  // Same RK4 momentum bound as R3; includes constant body-frame rotor momentum.
  EXPECT_LE((angular_momentum_N(x, c.mp) - H0).norm() / H0.norm(), 1e-10);
  // Omitting h_rotor from H_N must be distinguishable on this trajectory.
  const Vec3 wrong0 = quat_rotate(rb::quat(c.x0), c.mp.inertia * rb::omega(c.x0));
  const Vec3 wrongT = quat_rotate(rb::quat(x), c.mp.inertia * rb::omega(x));
  EXPECT_GT((wrongT - wrong0).norm(), .1);  // Rotor contribution is O(.5) kg m^2/s.
}

TEST(Rotational, NormDriftWithoutRenormalization) {
  const auto c = case_R2();
  for (auto method : {Integrator::Euler, Integrator::RK4, Integrator::RKMK4}) {
    const auto x = integrate_case(c, method, .01, false).x_final;
    const double norm = std::abs(rb::quat(x).norm() - 1);
    const double E = std::abs(rotational_energy(x, c.mp) / rotational_energy(c.x0, c.mp) - 1);
    const double H = (angular_momentum_N(x, c.mp) - angular_momentum_N(c.x0, c.mp)).norm() /
                     angular_momentum_N(c.x0, c.mp).norm();
    std::printf("[norm/conservation] R2 %s projection=false norm=%.12g E=%.12g H=%.12g\n", to_string(method),
                norm, E, H);
    if (method == Integrator::RKMK4) EXPECT_LE(norm, 1e-13);  // Task 2.2 Lie norm bound.
    if (method == Integrator::Euler) EXPECT_GT(norm, .01);    // Positive Euler radial error after 1000 steps.
  }
}

TEST(Rotational, R4FlipOccursAndConserves) {
  const auto c = case_R4();
  Vec<13> x = c.x0;
  const double E0 = rotational_energy(x, c.mp);
  const Vec3 H0 = angular_momentum_N(x, c.mp);
  auto f = [&](double, const Vec<13>& y) { return rigid_body_deriv(y, c.mp, {}, c.g_N); };
  bool flip = false;
  double maxE = 0, maxH = 0;
  StepStats stats;
  for (int i = 0; i < 60000; ++i) {
    x = step(Integrator::RK4, f, i * .001, x, .001, true, stats);
    flip = flip || rb::omega(x).y() < 0;
    maxE = std::max(maxE, std::abs(rotational_energy(x, c.mp) / E0 - 1));
    maxH = std::max(maxH, (angular_momentum_N(x, c.mp) - H0).norm() / H0.norm());
  }
  EXPECT_TRUE(flip);
  // Task 2.2: R4 demo preserves energy and momentum over the entire 60 s trajectory to 1e-8.
  EXPECT_LE(maxE, 1e-8);
  EXPECT_LE(maxH, 1e-8);
}
