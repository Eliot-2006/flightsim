#include "flightsim/aircraft/f16_model.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "flightsim/math/frames.hpp"
#include "flightsim/math/units.hpp"

namespace flightsim {
namespace {

double power_derivative(const F16Data& d, double power, double command) {
  double target = command;
  if (power >= d.power_threshold) {
    if (command < d.power_threshold) target = d.power_down_target;
    return d.power_high_rate * (target - power);
  }
  if (command >= d.power_threshold) target = d.power_up_target;
  const double delta = target - power;
  const double magnitude = std::abs(delta);
  const double rate = magnitude <= d.rtau_dp_low    ? d.rtau_low
                      : magnitude >= d.rtau_dp_high ? d.rtau_high
                                                    : d.rtau_mid_intercept + d.rtau_mid_slope * magnitude;
  return rate * delta;
}

double actuator_derivative(const F16Actuator& a, double position, double command) {
  return std::clamp(a.lag_rate * (std::clamp(command, a.pos_min, a.pos_max) - position), -a.rate_limit,
                    a.rate_limit);
}

}  // namespace

MassProps F16Model::mass_props() const {
  Mat3 inertia;
  inertia << data.Ixx, 0, -data.Ixz, 0, data.Iyy, 0, -data.Ixz, 0, data.Izz;
  return {data.mass, inertia, Vec3(include_rotor ? data.hx : 0, 0, 0)};
}

F16Model::Eval F16Model::evaluate(double t, const Vec<17>& x, const F16Input& u) const {
  if (!x.allFinite() || !std::isfinite(xcg) || !std::isfinite(u.throttle) || !std::isfinite(u.de) ||
      !std::isfinite(u.da) || !std::isfinite(u.dr)) {
    throw std::invalid_argument("F16Model: state, CG and inputs must be finite");
  }
  const auto& d = data;
  Eval e;
  e.surfaces << std::clamp(u.throttle, d.throttle_min, d.throttle_max), actuators ? x(ac::kDe) : u.de,
      actuators ? x(ac::kDa) : u.da, actuators ? x(ac::kDr) : u.dr;
  const double de = e.surfaces(1), da = e.surfaces(2) / d.aileron_norm, dr = e.surfaces(3) / d.rudder_norm;
  e.power_cmd = e.surfaces(0) <= d.gear_breakpoint
                    ? d.gear_low_slope * e.surfaces(0)
                    : d.gear_high_slope * e.surfaces(0) + d.gear_high_intercept;

  const Vec3 velocity_B = dcm_from_quat(rb::quat(x)).transpose() * (rb::vel(x) - wind.wind_N(t));
  const auto air = air_data_from_body_velocity(velocity_B);
  e.V = air.V;
  e.alpha = air.alpha;
  e.beta = air.beta;
  const double altitude = -x(rb::kPos + 2);
  const auto air_state = atmosphere(atm, altitude);
  e.mach = e.V / air_state.a;
  e.qbar = 0.5 * air_state.rho * e.V * e.V;

  // S&L uses rounded rtod for alpha/beta lookup degrees; the loader stores axes in exact SI radians.
  const double alpha = e.alpha * d.rtod * units::deg, beta = e.beta * d.rtod * units::deg;
  const double odd = beta > 0 ? 1.0 : beta < 0 ? -1.0 : 0.0;
  const Vec3 omega = rb::omega(x);
  const double inv_2V = e.V < 1e-9 ? 0 : 0.5 / e.V;
  const double pb = omega.x() * d.b * inv_2V, qc = omega.y() * d.cbar * inv_2V,
               rb_rate = omega.z() * d.b * inv_2V;

  const double cx = d.cx(alpha, de, monitor) + d.cxq(alpha, monitor) * qc;
  const double cy = d.cy_beta_gain * beta + d.cy_aileron_gain * da + d.cy_rudder_gain * dr +
                    d.cyr(alpha, monitor) * rb_rate + d.cyp(alpha, monitor) * pb;
  const double cz = d.cz(alpha, monitor) * (1 - std::pow(beta / d.cz_beta_scale, 2)) +
                    d.cz_elevator_gain * de / d.cz_elevator_scale + d.czq(alpha, monitor) * qc;
  const double cl = odd * d.cl(alpha, std::abs(beta), monitor) + d.dlda(alpha, beta, monitor) * da +
                    d.dldr(alpha, beta, monitor) * dr + d.clr(alpha, monitor) * rb_rate +
                    d.clp(alpha, monitor) * pb;
  const double cm = d.cm(alpha, de, monitor) + d.cmq(alpha, monitor) * qc + cz * (d.xcg_ref - xcg);
  const double cn = odd * d.cn(alpha, std::abs(beta), monitor) + d.dnda(alpha, beta, monitor) * da +
                    d.dndr(alpha, beta, monitor) * dr + d.cnr(alpha, monitor) * rb_rate +
                    d.cnp(alpha, monitor) * pb - cy * (d.xcg_ref - xcg) * d.cbar / d.b;

  const double military = d.thrust_mil(altitude, e.mach, monitor);
  const double power = x(ac::kPower);
  if (power < d.power_levels[1]) {
    const double idle = d.thrust_idle(altitude, e.mach, monitor);
    e.thrust =
        idle + (military - idle) * (power - d.power_levels[0]) / (d.power_levels[1] - d.power_levels[0]);
  } else {
    const double maximum = d.thrust_max(altitude, e.mach, monitor);
    e.thrust = military +
               (maximum - military) * (power - d.power_levels[1]) / (d.power_levels[2] - d.power_levels[1]);
  }
  const double qs = e.qbar * d.S;
  e.fm.F_B = qs * Vec3(cx, cy, cz) + Vec3(e.thrust, 0, 0);
  e.fm.M_B = qs * Vec3(d.b * cl, d.cbar * cm, d.b * cn);
  return e;
}

Vec<17> F16Model::deriv(double t, const Vec<17>& x, const F16Input& u) const {
  const auto e = evaluate(t, x, u);
  Vec<17> dx = Vec<17>::Zero();
  dx.head<13>() =
      rigid_body_deriv(x.head<13>(), mass_props(), e.fm, Vec3(0, 0, grav.gravity(-x(rb::kPos + 2))));
  dx(ac::kPower) = power_derivative(data, x(ac::kPower), e.power_cmd);
  if (actuators) {
    dx(ac::kDe) = actuator_derivative(data.elevator, x(ac::kDe), u.de);
    dx(ac::kDa) = actuator_derivative(data.aileron, x(ac::kDa), u.da);
    dx(ac::kDr) = actuator_derivative(data.rudder, x(ac::kDr), u.dr);
  }
  return dx;
}

}  // namespace flightsim
