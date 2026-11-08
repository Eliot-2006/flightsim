#pragma once

#include "flightsim/math/quat.hpp"

// Closed-form translational solutions (PRD 8.2). All in NED SI units, z down.
namespace flightsim::validation {

struct PosVel {
  Vec3 p, v;
};

// T1: constant acceleration a: p = p0 + v0 t + a t^2/2.
PosVel exact_constant_accel(const Vec3& p0, const Vec3& v0, const Vec3& a, double t);

// T2: m vdot = -c v + m g. v = v_inf + (v0 - v_inf) e^{-ct/m}, v_inf = m g / c; p integrated analytically.
PosVel exact_linear_drag(double m, double c, const Vec3& g, const Vec3& p0, const Vec3& v0, double t);

// T3: vertical fall from rest with F = -k |v| v, v_t = sqrt(m g/k):
// v = v_t tanh(g t/v_t), z = (v_t^2/g) ln cosh(g t/v_t) (down positive). Returned along the NED z axis.
PosVel exact_quadratic_drag_drop(double m, double k, double g, double t);

// T4: m pddot = -diag(k) p, per-axis w_i = sqrt(k_i/m).
PosVel exact_spring(double m, const Vec3& k_diag, const Vec3& p0, const Vec3& v0, double t);

// T5: body spins at constant rate Omega about its own z axis (torque-free isotropic body) with constant
// body-x force F. C_NB(t) = C_NB(0) Rz(Omega t), so with F_N(t) = C_NB(0) Rz(Omega t) [F,0,0]:
//   v = v0 + C_NB(0) (F/(m Omega)) [sin Wt, 1 - cos Wt, 0]
//   p = p0 + v0 t + C_NB(0) (F/(m Omega^2)) [1 - cos Wt, Wt - sin Wt, 0]
// q(t) = q0 (x) [cos(Wt/2), 0, 0, sin(Wt/2)]. For q0 = identity this is the plan's closed form.
struct RotatingForceSolution {
  PosVel pv;
  Vec4 q;
};
RotatingForceSolution exact_rotating_body_force(double m, double F, double Omega, const Vec4& q0,
                                                const Vec3& p0, const Vec3& v0, double t);

}  // namespace flightsim::validation
