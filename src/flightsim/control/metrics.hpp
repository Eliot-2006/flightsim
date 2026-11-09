#pragma once
#include <vector>

namespace flightsim {

// Control-performance metrics on sampled time histories (t strictly increasing). Pure functions: the runner
// feeds them its recorded samples; times are reported at sample instants (resolution = the sample interval).

struct StepMetrics {
  double step;             // y_final - y_initial (signed)
  double overshoot_pct;    // max excursion beyond y_final in the step direction, % of |step|; 0 if none
  double settling_2pct_s;  // from t_cmd to the first sample after which y stays within 2 % of |step|; NaN if
                           // never
  double steady_state_error;  // max |y - y_final| over the last `window_s` of the record
};
// Response y(t) to a command step at t_cmd from y_initial to y_final. Throws std::invalid_argument for a zero
// step, mismatched/short series or t_cmd outside the record.
StepMetrics step_metrics(const std::vector<double>& t, const std::vector<double>& y, double t_cmd,
                         double y_initial, double y_final, double window_s = 10.0);

struct DisturbanceMetrics {
  double peak_deviation;   // max |y - y_ref| over [t_from, t_to]
  double recovery_2pct_s;  // from t_from to the first sample after which |y - y_ref| stays <= 2 % of the peak
                           // until t_to; NaN if it has not recovered by t_to; 0 for a zero disturbance
};
DisturbanceMetrics disturbance_metrics(const std::vector<double>& t, const std::vector<double>& y,
                                       double y_ref, double t_from, double t_to);

// max |y - y_ref| over samples with t_from <= t <= t_to.
double peak_abs_deviation(const std::vector<double>& t, const std::vector<double>& y, double y_ref,
                          double t_from, double t_to);
// max |dy/dt| by first differences of consecutive samples (actuator rate).
double peak_abs_rate(const std::vector<double>& t, const std::vector<double>& y);

}  // namespace flightsim
