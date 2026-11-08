#pragma once

#include "flightsim/math/quat.hpp"

namespace flightsim {

struct AirData {
  double V, alpha, beta;
};

// v_air_B = (u, v, w) in body axes. V = |v_air|, alpha = atan2(w, u), beta = asin(v/V) (PRD 4.1).
// V < 1e-9 m/s: angles are undefined, alpha = beta = 0.
AirData air_data_from_body_velocity(const Vec3& v_air_B);

// [phi_dot, theta_dot, psi_dot] = E * omega_B for 3-2-1 Euler angles.
// Singular at theta = +-pi/2 (entries are non-finite there); callers must stay away from the pole.
Mat3 euler_rate_matrix(double phi, double theta);

}  // namespace flightsim
