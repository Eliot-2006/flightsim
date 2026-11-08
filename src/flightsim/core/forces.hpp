#pragma once

#include <functional>

#include "flightsim/core/rigid_body.hpp"

namespace flightsim {

using ForceModel = std::function<ForceMoment(double t, const Vec<13>& x)>;

ForceModel linear_drag(double c);       // F_N = -c v_N
ForceModel quadratic_drag(double k);    // F_N = -k |v_N| v_N
ForceModel spring(const Vec3& k_diag);  // F_N = -diag(k) p_N
ForceModel body_force(const Vec3& F_B);
ForceModel body_torque(const Vec3& M_B);

}  // namespace flightsim
