#include "flightsim/validation/cases.hpp"

#include "flightsim/math/units.hpp"
#include "flightsim/validation/analytic.hpp"

namespace flightsim::validation {

namespace {

Vec<13> make_state(const Vec3& p, const Vec3& v, const Vec4& q, const Vec3& w) {
  Vec<13> x;
  x << p, v, q, w;
  return x;
}

const Vec4 kIdentityQ(1, 0, 0, 0);

Vec<13> from_pv(const PosVel& s, const Vec4& q = kIdentityQ, const Vec3& w = Vec3::Zero()) {
  return make_state(s.p, s.v, q, w);
}

}  // namespace

// T0: m = 1 kg, no force, g = 0, v0 = (3, -2, 1) m/s, T = 5 s. p = p0 + v0 t.
CaseDef case_T0() {
  const Vec3 p0(1, 2, 3), v0(3, -2, 1);
  CaseDef c{"T0",
            MassProps{1.0, Mat3::Identity()},
            Vec3::Zero(),
            {},
            make_state(p0, v0, kIdentityQ, Vec3::Zero()),
            5.0,
            {}};
  c.exact = [=](double t) { return from_pv(exact_constant_accel(p0, v0, Vec3::Zero(), t)); };
  return c;
}

// T1: m = 1 kg, no force, uniform gravity g_N = (0, 0, g0), v0 = (1, 2, 3) m/s, T = 5 s.
CaseDef case_T1() {
  const Vec3 p0(0, 0, -100), v0(1, 2, 3), g(0, 0, units::g0);
  CaseDef c{"T1", MassProps{1.0, Mat3::Identity()}, g, {}, make_state(p0, v0, kIdentityQ, Vec3::Zero()), 5.0,
            {}};
  c.exact = [=](double t) { return from_pv(exact_constant_accel(p0, v0, g, t)); };
  return c;
}

// T2: m = 2 kg, linear drag c = 0.5 N s/m, v0 = (10, -5, -20) m/s, g = 9.80665 m/s^2 down, T = 5 s.
CaseDef case_T2() {
  const double m = 2.0, k = 0.5;
  const Vec3 p0 = Vec3::Zero(), v0(10, -5, -20), g(0, 0, units::g0);
  CaseDef c{"T2",
            MassProps{m, Mat3::Identity()},
            g,
            {linear_drag(k)},
            make_state(p0, v0, kIdentityQ, Vec3::Zero()),
            5.0,
            {}};
  c.exact = [=](double t) { return from_pv(exact_linear_drag(m, k, g, p0, v0, t)); };
  return c;
}

// T3: m = 1 kg, quadratic drag k = 0.01 kg/m (v_t = sqrt(m g/k) = 31.3 m/s), dropped from rest, T = 10 s.
CaseDef case_T3() {
  const double m = 1.0, k = 0.01;
  CaseDef c{"T3",
            MassProps{m, Mat3::Identity()},
            Vec3(0, 0, units::g0),
            {quadratic_drag(k)},
            make_state(Vec3::Zero(), Vec3::Zero(), kIdentityQ, Vec3::Zero()),
            10.0,
            {}};
  c.exact = [=](double t) { return from_pv(exact_quadratic_drag_drop(m, k, units::g0, t)); };
  return c;
}

// T4: m = 1 kg, spring k = (1, 4, 9) N/m, p0 = (1, 1, 1) m, v0 = (0, 0.5, 0) m/s, g = 0, T = 20 s.
CaseDef case_T4() {
  const double m = 1.0;
  const Vec3 k(1, 4, 9), p0(1, 1, 1), v0(0, 0.5, 0);
  CaseDef c{"T4",
            MassProps{m, Mat3::Identity()},
            Vec3::Zero(),
            {spring(k)},
            make_state(p0, v0, kIdentityQ, Vec3::Zero()),
            20.0,
            {}};
  c.exact = [=](double t) { return from_pv(exact_spring(m, k, p0, v0, t)); };
  return c;
}

// T5: m = 1 kg, I = I3, omega_B = (0, 0, 2) rad/s (torque-free about a principal axis of an isotropic body,
// so constant), body-x force F_B = (1, 0, 0) N, g = 0, T = 10 s. q0 is a generic tilted attitude (phi, theta,
// psi) = (0.3, -0.4, 0.5) rad so that C_NB vs C_BN *and* q (x) w vs w (x) q errors are both visible; for q0 =
// identity the closed form reduces to v_N = v0 + (F/(m Omega)) [sin Wt, 1 - cos Wt, 0].
CaseDef case_T5() {
  const double m = 1.0, F = 1.0, Omega = 2.0;
  const Vec4 q0 = quat_from_euler(0.3, -0.4, 0.5);
  const Vec3 p0(1, -1, 2), v0(0.5, 0.25, -0.75);
  CaseDef c{"T5",
            MassProps{m, Mat3::Identity()},
            Vec3::Zero(),
            {body_force(Vec3(F, 0, 0))},
            make_state(p0, v0, q0, Vec3(0, 0, Omega)),
            10.0,
            {}};
  c.exact = [=](double t) {
    const auto s = exact_rotating_body_force(m, F, Omega, q0, p0, v0, t);
    return from_pv(s.pv, s.q, Vec3(0, 0, Omega));
  };
  return c;
}

}  // namespace flightsim::validation
