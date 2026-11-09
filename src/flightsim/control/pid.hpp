#pragma once

#include <optional>

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
  // Same law, but back-calculation tracks `applied`, the output actually reaching the plant (e.g. a measured
  // actuator position, which also reflects rate limits and lag), instead of the command saturation.
  double update(double setpoint, double measurement, double dt, double umin, double umax, double applied);
  // Bumpless: the next update outputs exactly clamp(u0). With ki == 0 there is no integrator to carry the
  // offset, so only that first output is u0; afterwards the plain P(D) law applies with no stored bias.
  void reset(double u0 = 0);
  bool antiwindup = true;

 private:
  double step(double setpoint, double measurement, double dt, double umin, double umax,
              std::optional<double> applied);
  PidGains g_;
  double integral_ = 0, derivative_ = 0, measurement_ = 0, initial_output_ = 0;
  bool measured_ = false, resetting_ = false;
};
}  // namespace flightsim
