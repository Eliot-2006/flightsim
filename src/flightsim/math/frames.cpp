#include "flightsim/math/frames.hpp"

#include <algorithm>
#include <cmath>

namespace flightsim {

AirData air_data_from_body_velocity(const Vec3& v) {
  const double V = v.norm();
  if (V < 1e-9) return {V, 0.0, 0.0};
  return {V, std::atan2(v.z(), v.x()), std::asin(std::clamp(v.y() / V, -1.0, 1.0))};
}

Mat3 euler_rate_matrix(double phi, double theta) {
  const double sp = std::sin(phi), cp = std::cos(phi), tt = std::tan(theta), ct = std::cos(theta);
  Mat3 E;
  E << 1, sp * tt, cp * tt, 0, cp, -sp, 0, sp / ct, cp / ct;
  return E;
}

}  // namespace flightsim
