#include "flightsim/control/metrics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace flightsim {
namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

void check(const std::vector<double>& t, const std::vector<double>& y) {
  if (t.size() != y.size() || t.size() < 2) {
    throw std::invalid_argument("metrics: need matching time/value series of at least 2 samples");
  }
}

// Time (relative to t_from) of the first sample in [i0, i1] after which |dev| stays <= band through i1;
// NaN if the last sample is still outside. `dev` is evaluated per sample index.
template <class Dev>
double settle_time(const std::vector<double>& t, std::size_t i0, std::size_t i1, double band, double t_from,
                   Dev&& dev) {
  std::size_t last_outside = i1 + 1;  // none
  for (std::size_t i = i0; i <= i1; ++i) {
    if (std::abs(dev(i)) > band) {
      last_outside = i;
    }
  }
  if (last_outside == i1 + 1) {
    return 0.0;
  }
  return last_outside == i1 ? kNaN : t[last_outside + 1] - t_from;
}

std::size_t first_at_or_after(const std::vector<double>& t, double time) {
  return static_cast<std::size_t>(std::lower_bound(t.begin(), t.end(), time) - t.begin());
}
}  // namespace

StepMetrics step_metrics(const std::vector<double>& t, const std::vector<double>& y, double t_cmd,
                         double y_initial, double y_final, double window_s) {
  check(t, y);
  const double step = y_final - y_initial;
  if (step == 0.0 || !std::isfinite(step)) {
    throw std::invalid_argument("step_metrics: command step must be finite and nonzero");
  }
  const std::size_t i0 = first_at_or_after(t, t_cmd);
  if (i0 >= t.size() - 1) {
    throw std::invalid_argument("step_metrics: command time outside the record");
  }
  const std::size_t last = t.size() - 1;
  StepMetrics m{step, 0.0, 0.0, 0.0};
  double worst = 0;
  for (std::size_t i = i0; i <= last; ++i) {
    worst = std::max(worst, (y[i] - y_final) * (step > 0 ? 1.0 : -1.0));
  }
  m.overshoot_pct = 100 * worst / std::abs(step);
  m.settling_2pct_s =
      settle_time(t, i0, last, 0.02 * std::abs(step), t_cmd, [&](std::size_t i) { return y[i] - y_final; });
  for (std::size_t i = first_at_or_after(t, t.back() - window_s); i <= last; ++i) {
    m.steady_state_error = std::max(m.steady_state_error, std::abs(y[i] - y_final));
  }
  return m;
}

double peak_abs_deviation(const std::vector<double>& t, const std::vector<double>& y, double y_ref,
                          double t_from, double t_to) {
  check(t, y);
  double peak = 0;
  for (std::size_t i = first_at_or_after(t, t_from); i < t.size() && t[i] <= t_to; ++i) {
    peak = std::max(peak, std::abs(y[i] - y_ref));
  }
  return peak;
}

DisturbanceMetrics disturbance_metrics(const std::vector<double>& t, const std::vector<double>& y,
                                       double y_ref, double t_from, double t_to) {
  check(t, y);
  const std::size_t i0 = first_at_or_after(t, t_from);
  std::size_t i1 = i0;
  while (i1 + 1 < t.size() && t[i1 + 1] <= t_to) {
    ++i1;
  }
  const double peak = peak_abs_deviation(t, y, y_ref, t_from, t_to);
  return {peak, peak == 0.0 ? 0.0 : settle_time(t, i0, i1, 0.02 * peak, t_from, [&](std::size_t i) {
            return y[i] - y_ref;
          })};
}

double peak_abs_rate(const std::vector<double>& t, const std::vector<double>& y) {
  check(t, y);
  double peak = 0;
  for (std::size_t i = 1; i < t.size(); ++i) {
    peak = std::max(peak, std::abs((y[i] - y[i - 1]) / (t[i] - t[i - 1])));
  }
  return peak;
}
}  // namespace flightsim
