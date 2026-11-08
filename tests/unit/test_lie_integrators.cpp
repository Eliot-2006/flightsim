#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "flightsim/core/integrators.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/validation/cases.hpp"

using namespace flightsim;

namespace {
Vec<13> initial() {
  Vec<13> x = Vec<13>::Zero();
  rb::quat(x) = quat_from_euler(0.3, -0.4, 0.5);
  rb::omega(x) = Vec3(0.3, 0.4, 1.3);
  return x;
}

Vec<13> free_deriv(double, const Vec<13>& x) {
  return rigid_body_deriv(x, MassProps{1, Vec3(1, 2, 3).asDiagonal()}, {}, Vec3::Zero());
}

// Deliberately wrong left-Jacobian sign; this must fail even though constant omega passes.
Vec<13> wrong_sign_step(const Vec<13>& x, double h) {
  auto stage = [&](const Vec<13>& v, const Vec3& theta) {
    Vec<13> y = v;
    rb::quat(y) = quat_mul(rb::quat(x), quat_from_rotvec(theta));
    const Vec3 w = rb::omega(y);
    return std::pair<Vec<13>, Vec3>{free_deriv(0, y),
                                    w - 0.5 * theta.cross(w) + theta.cross(theta.cross(w)) / 12.0};
  };
  const auto [k1, a1] = stage(x, Vec3::Zero());
  const auto [k2, a2] = stage(x + h / 2 * k1, h / 2 * a1);
  const auto [k3, a3] = stage(x + h / 2 * k2, h / 2 * a2);
  const auto [k4, a4] = stage(x + h * k3, h * a3);
  Vec<13> y = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  rb::quat(y) = quat_mul(rb::quat(x), quat_from_rotvec(h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)));
  return y;
}
}  // namespace

TEST(LieIntegrators, LieEulerPreservesNorm) {
  Vec<13> x = initial();
  StepStats stats;
  auto f = [](double, const Vec<13>& y) {
    Vec<13> d = Vec<13>::Zero();
    const Vec3 w = rb::omega(y);
    rb::quat(d) = 0.5 * quat_mul(rb::quat(y), Vec4(0, w.x(), w.y(), w.z()));
    return d;
  };
  for (int i = 0; i < 1000000; ++i) x = step(Integrator::LieEuler, f, 0, x, 0.001, false, stats);
  // Task 2.1 bound after one million exponential updates, without state projection.
  EXPECT_NEAR(rb::quat(x).norm(), 1.0, 1e-13);
  EXPECT_EQ(stats.nfev, 1000000);
}

TEST(LieIntegrators, Rkmk4ExactForConstantOmega) {
  Vec<13> x = initial();
  rb::omega(x) = Vec3(0, 0, 1.3);
  const Vec4 q0 = rb::quat(x);
  StepStats stats;
  auto f = [](double, const Vec<13>& y) {
    return rigid_body_deriv(y, MassProps{1, Mat3::Identity()}, {}, Vec3::Zero());
  };
  for (int i = 0; i < 1000; ++i) x = step(Integrator::RKMK4, f, i * 0.01, x, 0.01, false, stats);
  // PRD 10.3: RKMK4 exact constant-rate attitude to 1e-12 rad after 1000 steps.
  EXPECT_LE(attitude_error(quat_mul(q0, quat_from_rotvec(Vec3(0, 0, 13))), rb::quat(x)), 1e-12);
  EXPECT_EQ(stats.nfev, 4000);
}

TEST(LieIntegrators, Rkmk4SignRegression) {
  const Vec<13> x = initial();
  Vec<13> ref = x;
  const double h = 0.01;
  for (int i = 0; i < 1000; ++i) ref = rk4_step(free_deriv, i * h / 1000, ref, h / 1000);
  StepStats stats;
  const Vec<13> good = step(Integrator::RKMK4, free_deriv, 0, x, h, false, stats);
  const double error = attitude_error(rb::quat(ref), rb::quat(good));
  const double wrong = attitude_error(rb::quat(ref), rb::quat(wrong_sign_step(x, h)));
  // Task 2.1: local O(h^5) error <= 1e-9 at h=.01; wrong sign >=100 times worse.
  EXPECT_LE(error, 1e-9);
  EXPECT_GE(wrong, 100 * error);
  EXPECT_GT(wrong, 1e-8);  // O(h^3) commutator error, far above reference roundoff.
}

TEST(LieIntegrators, SymplecticEulerOscillatorEnergyBounded) {
  const auto c = validation::case_T4();
  auto f = [&](double, const Vec<13>& x) { return rigid_body_deriv(x, c.mp, c.forces.front()(0, x), c.g_N); };
  auto energy = [](const Vec<13>& x) {
    return 0.5 * (rb::vel(x).squaredNorm() + Vec3(1, 4, 9).dot(rb::pos(x).cwiseAbs2()));
  };
  Vec<13> x = c.x0;
  const double e0 = energy(x);
  double max_drift = 0, first = 0, last = 0;
  StepStats stats;
  for (int i = 0; i < 100000; ++i) {
    x = step(Integrator::SymplecticEuler, f, i * .01, x, .01, true, stats);
    const double drift = (energy(x) - e0) / e0;
    max_drift = std::max(max_drift, std::abs(drift));
    if (i < 10000) first += drift / 10000;
    if (i >= 90000) last += drift / 10000;
  }
  // Task 2.1 energy envelope and first/last 10% means, at h=.01 for 1e5 steps.
  EXPECT_LE(max_drift, 0.05);
  EXPECT_LE(std::abs(last - first), 1e-3);
  EXPECT_EQ(stats.nfev, 100000);
}

TEST(LieIntegrators, StageAttitudeCouplesToVectorStatesAndExtraState) {
  Vec<14> x = Vec<14>::Zero();
  x.head<13>() = initial();
  x(13) = 2;
  auto f = [](double t, const Vec<14>& y) {
    Vec<14> d = Vec<14>::Zero();
    d.head<13>() = free_deriv(t, y.head<13>());
    rb::vel(d) = quat_rotate(rb::quat(y), Vec3::UnitX());
    d(13) = -y(13);
    rb::quat(d).setConstant(999);  // Lie stages must ignore the supplied quaternion derivative.
    return d;
  };
  auto reference_f = [&](double t, const Vec<14>& y) {
    Vec<14> d = f(t, y);
    rb::quat(d) = rb::quat(free_deriv(t, y.head<13>()));
    return d;
  };
  Vec<14> ref = x;
  for (int i = 0; i < 1000; ++i) ref = rk4_step(reference_f, i * 1e-5, ref, 1e-5);
  StepStats stats;
  const Vec<14> y = step(Integrator::RKMK4, f, 0, x, .01, false, stats);
  // Local RK4 O(h^5) at h=.01; 1e-9 is the sign-regression local-error bound.
  EXPECT_LE((y - ref).norm(), 1e-9);
  EXPECT_EQ(stats.nfev, 4);
}
