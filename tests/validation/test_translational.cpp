#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/validation/analytic.hpp"
#include "flightsim/validation/cases.hpp"
#include "flightsim/validation/convergence.hpp"

using namespace flightsim;
using namespace flightsim::validation;

namespace {

// Observed order of one error component over a PRD 6.3 study (h0 = 0.1, 8 levels).
struct Orders {
  OrderFit pos, vel;
};

Orders study_orders(const CaseDef& c, Integrator integ) {
  const auto rows = convergence_study(c, integ, 0.1, 8);
  std::vector<double> h, ep, ev;
  for (const auto& r : rows) {
    h.push_back(r.h);
    ep.push_back(r.err.pos);
    ev.push_back(r.err.vel);
  }
  const Vec<13> xT = c.exact(c.t_final);
  return {observed_order(h, ep, xT.segment<3>(rb::kPos).norm()),
          observed_order(h, ev, xT.segment<3>(rb::kVel).norm())};
}

double spring_energy(const Vec<13>& x, double m, const Vec3& k) {
  return 0.5 * m * rb::vel(x).squaredNorm() + 0.5 * k.dot(rb::pos(x).cwiseProduct(rb::pos(x)));
}

}  // namespace

TEST(Translational, T0ZeroForceExactAllIntegrators) {
  const CaseDef c = case_T0();
  for (Integrator i : {Integrator::Euler, Integrator::RK4}) {
    const RunResult r = integrate_case(c, i, 0.05);
    const ErrorNorms e = state_error(r.x_final, c.exact(c.t_final));
    // Tolerance 1e-12 (plan value): p + n*h*v accumulates n=~100 roundoffs of ~1e-16.
    EXPECT_LE(e.pos, 1e-12) << to_string(i);
    EXPECT_LE(e.vel, 1e-12) << to_string(i);
    EXPECT_LE(e.att, 1e-12) << to_string(i);
  }
}

TEST(Translational, T1ConstantAccelRk4Exact) {
  const CaseDef c = case_T1();
  const RunResult r = integrate_case(c, Integrator::RK4, 0.05);
  const ErrorNorms e = state_error(r.x_final, c.exact(c.t_final));
  // Tolerance 1e-10 (plan value): RK4 integrates a quadratic exactly; only roundoff over ~100 steps remains.
  EXPECT_LE(e.pos, 1e-10);
  EXPECT_LE(e.vel, 1e-10);
}

// PRD 10.3 / plan: first-order methods in [0.95, 1.05], fourth-order in [3.9, 4.1].
TEST(Translational, T2EulerOrder1) {
  const Orders o = study_orders(case_T2(), Integrator::Euler);
  EXPECT_GE(o.pos.order, 0.95);
  EXPECT_LE(o.pos.order, 1.05);
  EXPECT_GE(o.vel.order, 0.95);
  EXPECT_LE(o.vel.order, 1.05);
  std::printf("[order] T2 Euler pos=%.4f vel=%.4f (levels %d-%d / %d-%d)\n", o.pos.order, o.vel.order,
              o.pos.first_level, o.pos.last_level, o.vel.first_level, o.vel.last_level);
}

TEST(Translational, T2Rk4Order4) {
  const Orders o = study_orders(case_T2(), Integrator::RK4);
  EXPECT_GE(o.pos.order, 3.9);
  EXPECT_LE(o.pos.order, 4.1);
  EXPECT_GE(o.vel.order, 3.9);
  EXPECT_LE(o.vel.order, 4.1);
  std::printf("[order] T2 RK4 pos=%.4f vel=%.4f (levels %d-%d / %d-%d)\n", o.pos.order, o.vel.order,
              o.pos.first_level, o.pos.last_level, o.vel.first_level, o.vel.last_level);
}

TEST(Translational, T3Rk4Order4) {
  const Orders o = study_orders(case_T3(), Integrator::RK4);
  EXPECT_GE(o.pos.order, 3.9);
  EXPECT_LE(o.pos.order, 4.1);
  EXPECT_GE(o.vel.order, 3.9);
  EXPECT_LE(o.vel.order, 4.1);
  std::printf("[order] T3 RK4 pos=%.4f vel=%.4f (levels %d-%d / %d-%d)\n", o.pos.order, o.vel.order,
              o.pos.first_level, o.pos.last_level, o.vel.first_level, o.vel.last_level);
}

TEST(Translational, T3EulerOrder1) {
  const Orders o = study_orders(case_T3(), Integrator::Euler);
  EXPECT_GE(o.pos.order, 0.95);
  EXPECT_LE(o.pos.order, 1.05);
  EXPECT_GE(o.vel.order, 0.95);
  EXPECT_LE(o.vel.order, 1.05);
  std::printf("[order] T3 Euler pos=%.4f vel=%.4f (levels %d-%d / %d-%d)\n", o.pos.order, o.vel.order,
              o.pos.first_level, o.pos.last_level, o.vel.first_level, o.vel.last_level);
}

TEST(Translational, T4EulerEnergyGrows) {
  const CaseDef c = case_T4();
  const Vec3 k(1, 4, 9);
  const double h = 0.01;
  const Vec<13> x0 = c.x0;
  const RunResult r = integrate_case(c, Integrator::Euler, h);
  const double E0 = spring_energy(x0, 1.0, k), ET = spring_energy(r.x_final, 1.0, k);
  EXPECT_GT(ET, E0);
  // Explicit Euler multiplies each mode's energy by exactly (1 + w^2 h^2) per step (algebraic identity),
  // modes decouple.
  const long n = std::lround(c.t_final / h);
  const Vec3 p0 = rb::pos(x0), v0 = rb::vel(x0);
  double expected = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double w2 = k(i);
    expected +=
        (0.5 * v0(i) * v0(i) + 0.5 * w2 * p0(i) * p0(i)) * std::pow(1.0 + w2 * h * h, static_cast<double>(n));
  }
  // Tolerance 1e-9 rel: n = 2000 steps of roundoff (~2000 * 1e-16 = 2e-13) times a 5000x safety for pow()
  // rounding.
  EXPECT_NEAR(ET / expected, 1.0, 1e-9);
  // Slowest mode alone (w = 1) grows by (1+1e-4)^2000 = 1.22; the total must exceed that.
  EXPECT_GT(ET / E0, 1.2);
}

TEST(Translational, T4Rk4EnergyDriftSmall) {
  const CaseDef c = case_T4();
  const Vec3 k(1, 4, 9);
  const RunResult r = integrate_case(c, Integrator::RK4, 0.01);
  const double E0 = spring_energy(c.x0, 1.0, k), ET = spring_energy(r.x_final, 1.0, k);
  const double drift = std::abs(ET - E0) / E0;
  // Plan bound 1e-6. Analytic estimate for context: RK4 damps a mode by ~(wh)^6/144 per step -> ~1e-8 over
  // 2000 steps.
  EXPECT_LE(drift, 1e-6);
  std::printf("[T4] RK4 h=0.01 relative energy drift = %.3e\n", drift);
}

namespace {

// Deliberately wrong derivative variants for T5, used to show what the test catches.
enum class Bug { None, TransposedDcm, LeftQuatMultiply };

Vec<13> t5_deriv(const CaseDef& c, Bug bug, const Vec<13>& x) {
  ForceMoment fm;
  for (const auto& f : c.forces) fm += f(0.0, x);
  Vec<13> d = rigid_body_deriv(x, c.mp, fm, c.g_N);
  if (bug == Bug::TransposedDcm) {
    rb::vel(d) = dcm_from_quat(rb::quat(x)).transpose() * fm.F_B / c.mp.mass;  // C_BN instead of C_NB
  } else if (bug == Bug::LeftQuatMultiply) {
    const Vec3 w = rb::omega(x);
    rb::quat(d) =
        0.5 * quat_mul(Vec4(0, w.x(), w.y(), w.z()), rb::quat(x));  // [0,w] (x) q instead of q (x) [0,w]
  }
  return d;
}

ErrorNorms run_t5(const CaseDef& c, Bug bug, double h) {
  Vec<13> x = c.x0;
  const long n = std::lround(c.t_final / h);
  for (long i = 0; i < n; ++i) {
    x = rk4_step<13>([&](double, const Vec<13>& y) { return t5_deriv(c, bug, y); }, 0.0, x, h);
    renormalize_quat(x);
  }
  return state_error(x, c.exact(c.t_final));
}

}  // namespace

TEST(Translational, T5DirectionOfRotation) {
  const CaseDef c = case_T5();
  // Correct model, RK4 h = 1e-3: error = O(T*(Omega*h)^4*Omega/120) ~ 1e-12; 1e-9 is the plan bound (~1000x
  // margin).
  const ErrorNorms good = run_t5(c, Bug::None, 1e-3);
  EXPECT_LE(good.pos, 1e-9);
  EXPECT_LE(good.vel, 1e-9);
  EXPECT_LE(good.att, 1e-9);
  // Same run through the public harness agrees (so integrate_case is the model the test exercises).
  const ErrorNorms via_harness =
      state_error(integrate_case(c, Integrator::RK4, 1e-3).x_final, c.exact(c.t_final));
  EXPECT_LE(via_harness.vel, 1e-9);

  // Transposed DCM (C_BN used where C_NB belongs): the force rotates the wrong way. Plan threshold 0.1 m/s.
  const ErrorNorms transposed = run_t5(c, Bug::TransposedDcm, 1e-3);
  EXPECT_GT(transposed.vel, 0.1);
  // Left-multiplied quaternion kinematics: the body spins about the wrong axis (q0 is tilted so this is
  // visible).
  const ErrorNorms left = run_t5(c, Bug::LeftQuatMultiply, 1e-3);
  EXPECT_GT(left.vel, 0.1);
  std::printf("[T5] correct vel err=%.2e; transposed-DCM vel err=%.3f; left-mult vel err=%.3f\n", good.vel,
              transposed.vel, left.vel);
}

TEST(Translational, ObservedOrderFitExcludesRoundoffFloor) {
  // Synthetic: e = 3 h^2 (exact power law, order 2) until it drops below the floor 1e3*eps*scale = 2.2e-13,
  // then a non-monotone roundoff plateau at ~1e-15. scale = 1. 26 levels reach the floor at k = 19.
  const int levels = 26;
  const double floor_ = 1e3 * 2.220446049250313e-16;
  std::vector<double> h, e;
  for (int k = 0; k < levels; ++k) {
    h.push_back(0.1 / std::pow(2.0, k));
    const double exact = 3.0 * h.back() * h.back();
    e.push_back(exact > floor_ ? exact : 1e-15 * (1.0 + 0.5 * (k % 2)));
  }
  const OrderFit fit = observed_order(h, e, 1.0);
  // Exact power law above the floor -> slope 2 up to fit roundoff (1e-9).
  EXPECT_NEAR(fit.order, 2.0, 1e-9);
  // The window's last level is the finest one still above the floor; no plateau level is included.
  EXPECT_GT(e[static_cast<std::size_t>(fit.last_level)], floor_);
  EXPECT_LT(fit.last_level, levels - 1);
  EXPECT_EQ(fit.first_level, 0);

  // Pre-asymptotic coarse levels are excluded by the 10% successive-order rule: corrupt the two coarsest
  // errors.
  std::vector<double> e2;
  for (std::size_t k = 0; k < 10; ++k) e2.push_back(3.0 * h[k] * h[k]);
  e2[0] *= 5.0;  // successive orders p0 differs wildly from p1
  const OrderFit fit2 = observed_order(std::vector<double>(h.begin(), h.begin() + 10), e2, 1.0);
  EXPECT_NEAR(fit2.order, 2.0, 1e-9);
  EXPECT_GE(fit2.first_level, 1);

  // Too few valid levels -> no fit reported.
  const OrderFit none = observed_order({0.1, 0.05}, {1e-20, 1e-21}, 1.0);
  EXPECT_TRUE(std::isnan(none.order));
}
