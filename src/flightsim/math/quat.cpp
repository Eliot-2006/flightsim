#include "flightsim/math/quat.hpp"

#include <algorithm>
#include <cmath>

namespace flightsim {

namespace {
constexpr double kSmallAngle = 1e-4;  // Taylor branches; truncation O(angle^6) < 1e-24 here
}

Vec4 quat_mul(const Vec4& a, const Vec4& b) {
  const Vec3 av = a.tail<3>(), bv = b.tail<3>();
  Vec4 r;
  r(0) = a(0) * b(0) - av.dot(bv);
  r.tail<3>() = a(0) * bv + b(0) * av + av.cross(bv);
  return r;
}

Vec4 quat_conj(const Vec4& q) { return Vec4(q(0), -q(1), -q(2), -q(3)); }

Vec4 quat_normalized(const Vec4& q) { return q / q.norm(); }

Mat3 dcm_from_quat(const Vec4& q) {
  const double q0 = q(0), q1 = q(1), q2 = q(2), q3 = q(3);
  Mat3 C;
  C << q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2),
      2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1),
      2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
  return C;
}

Vec3 quat_rotate(const Vec4& q, const Vec3& v_B) { return dcm_from_quat(q) * v_B; }

Vec4 quat_from_euler(double phi, double theta, double psi) {
  const Vec4 qz(std::cos(psi / 2), 0, 0, std::sin(psi / 2));
  const Vec4 qy(std::cos(theta / 2), 0, std::sin(theta / 2), 0);
  const Vec4 qx(std::cos(phi / 2), std::sin(phi / 2), 0, 0);
  return quat_mul(qz, quat_mul(qy, qx));
}

Vec3 euler_from_quat(const Vec4& q) {
  const Mat3 C = dcm_from_quat(q);
  const double cos_theta = std::hypot(C(2, 1), C(2, 2));
  const double theta =
      std::atan2(-C(2, 0), cos_theta);  // well-conditioned for all theta, no asin clamp needed
  if (cos_theta < 1e-12) return Vec3(0.0, theta, std::atan2(-C(0, 1), C(1, 1)));
  return Vec3(std::atan2(C(2, 1), C(2, 2)), theta, std::atan2(C(1, 0), C(0, 0)));
}

Vec4 quat_from_rotvec(const Vec3& rv) {
  const double a = rv.norm();
  double c, s_over_a;  // cos(a/2), sin(a/2)/a
  if (a < kSmallAngle) {
    const double a2 = a * a;
    c = 1.0 - a2 / 8.0 + a2 * a2 / 384.0;
    s_over_a = 0.5 - a2 / 48.0 + a2 * a2 / 3840.0;
  } else {
    c = std::cos(a / 2);
    s_over_a = std::sin(a / 2) / a;
  }
  return Vec4(c, s_over_a * rv.x(), s_over_a * rv.y(), s_over_a * rv.z());
}

Vec3 rotvec_from_quat(const Vec4& q_in) {
  const Vec4 q = q_in(0) < 0 ? Vec4(-q_in) : q_in;
  const Vec3 v = q.tail<3>();
  const double s = v.norm();
  double angle_over_s;  // 2*atan2(s, w)/s
  if (s < kSmallAngle) {
    const double r2 = (s / q(0)) * (s / q(0));
    angle_over_s = 2.0 / q(0) * (1.0 - r2 / 3.0 + r2 * r2 / 5.0);
  } else {
    angle_over_s = 2.0 * std::atan2(s, q(0)) / s;
  }
  return angle_over_s * v;
}

double attitude_error(const Vec4& q_ref, const Vec4& q) {
  const Vec4 d = quat_mul(quat_conj(q_ref), q);
  return 2.0 * std::atan2(d.tail<3>().norm(), std::abs(d(0)));
}

Mat3 skew(const Vec3& v) {
  Mat3 S;
  S << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
  return S;
}

}  // namespace flightsim
