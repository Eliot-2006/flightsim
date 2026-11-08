#include "flightsim/core/forces.hpp"

namespace flightsim {

ForceModel linear_drag(double c) {
  return [c](double, const Vec<13>& x) {
    ForceMoment fm;
    fm.F_N = -c * rb::vel(x);
    return fm;
  };
}

ForceModel quadratic_drag(double k) {
  return [k](double, const Vec<13>& x) {
    const Vec3 v = rb::vel(x);
    ForceMoment fm;
    fm.F_N = -k * v.norm() * v;
    return fm;
  };
}

ForceModel spring(const Vec3& k_diag) {
  return [k_diag](double, const Vec<13>& x) {
    ForceMoment fm;
    fm.F_N = -k_diag.cwiseProduct(rb::pos(x));
    return fm;
  };
}

ForceModel body_force(const Vec3& F_B) {
  return [F_B](double, const Vec<13>&) {
    ForceMoment fm;
    fm.F_B = F_B;
    return fm;
  };
}

ForceModel body_torque(const Vec3& M_B) {
  return [M_B](double, const Vec<13>&) {
    ForceMoment fm;
    fm.M_B = M_B;
    return fm;
  };
}

}  // namespace flightsim
