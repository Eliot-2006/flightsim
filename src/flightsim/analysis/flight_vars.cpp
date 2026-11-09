#include "flightsim/analysis/flight_vars.hpp"

#include <cmath>
#include <stdexcept>

#include "flightsim/math/frames.hpp"

namespace flightsim {

Vec<17> x17_from_x13(const F16Model& m, const Vec<13>& x13, const F16Input& u, double t) {
  Vec<17> x = Vec<17>::Zero();
  rb::pos(x) = Vec3(x13(9), x13(10), -x13(11));
  rb::quat(x) = quat_from_euler(x13(3), x13(4), x13(5));
  const Vec3 v_B = x13(0) * Vec3(std::cos(x13(1)) * std::cos(x13(2)), std::sin(x13(2)),
                                 std::sin(x13(1)) * std::cos(x13(2)));
  rb::vel(x) = quat_rotate(rb::quat(x), v_B) + m.wind.wind_N(t);
  rb::omega(x) = x13.segment<3>(6);
  x(ac::kPower) = x13(12);
  x.tail<3>() = Vec3(u.de, u.da, u.dr);
  return x;
}

Vec<13> x13_from_x17(const F16Model& m, const Vec<17>& x, double t) {
  const Vec3 v_B = dcm_from_quat(rb::quat(x)).transpose() * (rb::vel(x) - m.wind.wind_N(t));
  const auto air = air_data_from_body_velocity(v_B);
  Vec<13> x13;
  x13 << air.V, air.alpha, air.beta, euler_from_quat(rb::quat(x)), rb::omega(x), x(0), x(1), -x(2),
      x(ac::kPower);
  return x13;
}

Vec<13> f13(const F16Model& m, const Vec<13>& x13, const F16Input& u, double t) {
  const auto x = x17_from_x13(m, x13, u, t);
  const auto dx = m.deriv(t, x, u);
  const Mat3 C_BN = dcm_from_quat(rb::quat(x)).transpose();
  const Vec3 v_B = C_BN * (rb::vel(x) - m.wind.wind_N(t));
  const Vec3 dv_B = C_BN * rb::vel(dx) - rb::omega(x).cross(v_B);
  Vec<13> d = Vec<13>::Zero();
  const double V = v_B.norm();
  if (V >= 1e-9) {
    if (std::abs(std::cos(x13(2))) < 1e-9) {
      throw std::domain_error("f13: sideslip chart pole");
    }
    d(0) = v_B.dot(dv_B) / V;
    d(1) = (v_B.x() * dv_B.z() - v_B.z() * dv_B.x()) / (v_B.x() * v_B.x() + v_B.z() * v_B.z());
    d(2) = (V * dv_B.y() - v_B.y() * d(0)) / (V * V * std::cos(x13(2)));
  }
  if (std::abs(std::cos(x13(4))) < 1e-9) {
    throw std::domain_error("f13: Euler pitch chart pole");
  }
  d.segment<3>(3) = euler_rate_matrix(x13(3), x13(4)) * rb::omega(x);
  d.segment<3>(6) = rb::omega(dx);
  d.segment<3>(9) = Vec3(dx(0), dx(1), -dx(2));
  d(12) = dx(ac::kPower);
  return d;
}

}  // namespace flightsim
