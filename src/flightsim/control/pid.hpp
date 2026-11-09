#pragma once

namespace flightsim {
struct PidGains {
  double kp, ki, kd, tf, tt, umin, umax;
};
class Pid {
 public:
  explicit Pid(PidGains);
  double update(double setpoint, double measurement, double dt);
  // Same law with the output limits replaced for this sample, for a downstream stage tighter than the gains'
  // own limits; back-calculation then tracks that tighter bound.
  double update(double setpoint, double measurement, double dt, double umin, double umax);
  void reset(double u0 = 0);
  bool antiwindup = true;

 private:
  PidGains g_;
  double integral_ = 0, derivative_ = 0, measurement_ = 0, initial_output_ = 0;
  bool measured_ = false, resetting_ = false;
};
}  // namespace flightsim
