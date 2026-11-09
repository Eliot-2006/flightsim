#include "flightsim/validation/analytic.hpp"

#include <cmath>
#include <stdexcept>

#include "flightsim/math/elliptic.hpp"

namespace flightsim::validation {

PosVel exact_constant_accel(const Vec3& p0, const Vec3& v0, const Vec3& a, double t) {
  return {p0 + v0 * t + 0.5 * a * t * t, v0 + a * t};
}

PosVel exact_linear_drag(double m, double c, const Vec3& g, const Vec3& p0, const Vec3& v0, double t) {
  const Vec3 v_inf = m * g / c;
  const double e = std::exp(-c * t / m);
  return {p0 + v_inf * t + (v0 - v_inf) * (m / c) * (1.0 - e), v_inf + (v0 - v_inf) * e};
}

PosVel exact_quadratic_drag_drop(double m, double k, double g, double t) {
  const double vt = std::sqrt(m * g / k);
  const double a = g * t / vt;
  return {Vec3(0, 0, vt * vt / g * std::log(std::cosh(a))), Vec3(0, 0, vt * std::tanh(a))};
}

PosVel exact_spring(double m, const Vec3& k_diag, const Vec3& p0, const Vec3& v0, double t) {
  PosVel r;
  for (int i = 0; i < 3; ++i) {
    const double w = std::sqrt(k_diag(i) / m);
    r.p(i) = p0(i) * std::cos(w * t) + v0(i) / w * std::sin(w * t);
    r.v(i) = -p0(i) * w * std::sin(w * t) + v0(i) * std::cos(w * t);
  }
  return r;
}

RotatingForceSolution exact_rotating_body_force(double m, double F, double Omega, const Vec4& q0,
                                                const Vec3& p0, const Vec3& v0, double t) {
  const Mat3 C0 = dcm_from_quat(q0);
  const double s = std::sin(Omega * t), c = std::cos(Omega * t);
  const Vec3 dv = F / (m * Omega) * Vec3(s, 1.0 - c, 0.0);
  const Vec3 dp = F / (m * Omega * Omega) * Vec3(1.0 - c, Omega * t - s, 0.0);
  return {{p0 + v0 * t + C0 * dp, v0 + C0 * dv},
          quat_mul(q0, Vec4(std::cos(0.5 * Omega * t), 0.0, 0.0, std::sin(0.5 * Omega * t)))};
}

double rotational_energy(const Vec<13>& x, const MassProps& mp) {
  return 0.5 * rb::omega(x).dot(mp.inertia * rb::omega(x));
}

Vec3 angular_momentum_N(const Vec<13>& x, const MassProps& mp) {
  return quat_rotate(rb::quat(x), mp.inertia * rb::omega(x) + mp.h_rotor);
}

Vec<13> exact_symmetric_top(const Vec<13>& x0, const MassProps& mp, double t) {
  const Vec3 w0 = rb::omega(x0);
  const double A = mp.inertia(0, 0), C = mp.inertia(2, 2);
  const Vec3 H = angular_momentum_N(x0, mp);
  const double lambda = w0.z() * (C - A) / A;
  Vec<13> x = x0;
  rb::pos(x) += t * rb::vel(x0);
  rb::quat(x) = quat_mul(quat_mul(quat_from_rotvec(t * H / A), rb::quat(x0)),
                         quat_from_rotvec(Vec3(0, 0, -lambda * t)));
  rb::omega(x) = Vec3(w0.x() * std::cos(lambda * t) - w0.y() * std::sin(lambda * t),
                      w0.x() * std::sin(lambda * t) + w0.y() * std::cos(lambda * t), w0.z());
  return x;
}

Vec3 exact_asymmetric_omega(const Vec3& I, const Vec3& w0, double t) {
  const double A = I.x(), B = I.y(), C = I.z();
  const double twiceE = w0.dot(I.cwiseProduct(w0)), L2 = I.cwiseProduct(w0).squaredNorm();
  if (!(0 < A && A < B && B < C && L2 > twiceE * B) || w0.y() != 0 || !std::isfinite(t))
    throw std::invalid_argument("LL reference requires 0<I1<I2<I3, L^2>2E I2, omega2(0)=0, finite t");
  const double frequency = std::sqrt((C - B) * (L2 - twiceE * A) / (A * B * C));
  const double m = (B - A) * (twiceE * C - L2) / ((C - B) * (L2 - twiceE * A));
  const auto [sn, cn, dn] = jacobi_elliptic(frequency * t, m);
  const double a2 = std::copysign(std::sqrt((twiceE * C - L2) / (B * (C - B))), w0.x() * w0.z());
  // omega2 has the sign of omega1(0)*omega3(0): Euler gives omega2dot(0)=(C-A)w1w3/B.
  return Vec3(w0.x() * cn, a2 * sn, w0.z() * dn);
}
}  // namespace flightsim::validation
