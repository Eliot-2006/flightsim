#include "flightsim/control/pid.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace flightsim {
Pid::Pid(PidGains g) : g_(g) {
  for (double v : {g.kp, g.ki, g.kd, g.tf, g.tt, g.umin, g.umax}) {
    if (!std::isfinite(v)) {
      throw std::invalid_argument("PID gains and limits must be finite");
    }
  }
  if (g.tf < 0 || g.tt <= 0 || g.umin > g.umax) {
    throw std::invalid_argument("PID requires tf>=0, tt>0 and umin<=umax");
  }
}
double Pid::update(double setpoint, double measurement, double dt) {
  return update(setpoint, measurement, dt, g_.umin, g_.umax);
}
double Pid::update(double setpoint, double measurement, double dt, double umin, double umax) {
  if (!std::isfinite(setpoint) || !std::isfinite(measurement) || !std::isfinite(dt) || dt <= 0 ||
      !std::isfinite(umin) || !std::isfinite(umax) || umin > umax) {
    throw std::invalid_argument("PID requires finite inputs, dt>0 and umin<=umax");
  }
  const double error = setpoint - measurement;
  const double proportional = g_.kp * error;
  if (resetting_) {
    integral_ = initial_output_ - proportional;
    derivative_ = 0;
  } else {
    derivative_ = measured_ ? (g_.tf * derivative_ + measurement - measurement_) / (g_.tf + dt) : 0;
    integral_ += dt * g_.ki * error;
  }
  measurement_ = measurement;
  measured_ = true;
  const double raw = proportional + integral_ - g_.kd * derivative_;
  // Return the requested value directly: subtracting/re-adding P can lose low output bits.
  const double output = resetting_ ? initial_output_ : std::clamp(raw, umin, umax);
  resetting_ = false;
  if (antiwindup) {
    integral_ += dt * (output - raw) / g_.tt;
  }
  return output;
}
void Pid::reset(double u0) {
  if (!std::isfinite(u0)) {
    throw std::invalid_argument("PID reset output must be finite");
  }
  initial_output_ = std::clamp(u0, g_.umin, g_.umax);
  resetting_ = true;
  measured_ = false;
}
}  // namespace flightsim
