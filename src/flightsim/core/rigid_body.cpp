#include "flightsim/core/rigid_body.hpp"

namespace flightsim {

Mat3 MassProps::inertia_inv() const { return inertia.inverse(); }

ForceMoment& operator+=(ForceMoment& a, const ForceMoment& b) {
  a.F_N += b.F_N;
  a.F_B += b.F_B;
  a.M_B += b.M_B;
  return a;
}

Vec<13> rigid_body_deriv(const Vec<13>& x, const MassProps& mp, const ForceMoment& fm, const Vec3& g_N) {
  const Vec4 q = rb::quat(x);
  const Vec3 w = rb::omega(x);
  Vec<13> d;
  rb::pos(d) = rb::vel(x);
  rb::vel(d) = (fm.F_N + quat_rotate(q, fm.F_B)) / mp.mass + g_N;
  rb::quat(d) = 0.5 * quat_mul(q, Vec4(0.0, w.x(), w.y(), w.z()));
  rb::omega(d) = mp.inertia_inv() * (fm.M_B - w.cross(mp.inertia * w + mp.h_rotor));
  return d;
}

}  // namespace flightsim
