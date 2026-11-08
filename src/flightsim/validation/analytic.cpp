#include "flightsim/validation/analytic.hpp"

#include <cmath>

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

}  // namespace flightsim::validation
