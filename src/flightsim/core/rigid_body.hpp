#pragma once

#include "flightsim/core/state.hpp"
#include "flightsim/math/quat.hpp"

namespace flightsim {

struct MassProps {
  double mass;
  Mat3 inertia;                 // about the CG, body axes; products of inertia allowed
  Vec3 h_rotor = Vec3::Zero();  // constant rotor angular momentum, body axes
  Mat3 inertia_inv() const;
};

// Non-gravitational force/moment. F_N is expressed in NED, F_B and M_B in body axes (moment about the CG).
struct ForceMoment {
  Vec3 F_N = Vec3::Zero();
  Vec3 F_B = Vec3::Zero();
  Vec3 M_B = Vec3::Zero();
};
ForceMoment& operator+=(ForceMoment& a, const ForceMoment& b);

// PRD 5.1. g_N = gravity vector in NED. Pure function of its arguments.
//   pdot = v;  vdot = (F_N + C_NB F_B)/m + g_N;  qdot = 1/2 q (x) [0, w_B];
//   wdot = I^-1 (M_B - w x (I w + h_rotor))
Vec<13> rigid_body_deriv(const Vec<13>& x, const MassProps& mp, const ForceMoment& fm, const Vec3& g_N);

}  // namespace flightsim
