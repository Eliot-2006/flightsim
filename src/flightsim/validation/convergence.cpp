#include "flightsim/validation/convergence.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "flightsim/math/quat.hpp"

namespace flightsim::validation {

ErrorNorms state_error(const Vec<13>& x, const Vec<13>& ref) {
  return {(rb::pos(x) - rb::pos(ref)).norm(), (rb::vel(x) - rb::vel(ref)).norm(),
          attitude_error(Vec4(rb::quat(ref)), Vec4(rb::quat(x))), (rb::omega(x) - rb::omega(ref)).norm()};
}

RunResult integrate_case(const CaseDef& c, Integrator integ, double h, bool renormalize) {
  const long n = std::lround(c.t_final / h);
  if (n < 1 || std::abs(static_cast<double>(n) * h - c.t_final) > 1e-9 * c.t_final)
    throw std::invalid_argument(c.name + ": step h does not divide t_final");
  auto f = [&c](double t, const Vec<13>& x) {
    ForceMoment fm;
    for (const auto& force : c.forces) fm += force(t, x);
    return rigid_body_deriv(x, c.mp, fm, c.g_N);
  };
  Vec<13> x = c.x0;
  StepStats stats;
  for (long i = 0; i < n; ++i) x = step(integ, f, static_cast<double>(i) * h, x, h, renormalize, stats);
  return {x, stats.nfev};
}

std::vector<ConvergenceRow> convergence_study(const CaseDef& c, Integrator integ, double h0, int levels) {
  if (!c.exact)
    throw std::invalid_argument(c.name + ": no closed form; fine-reference studies arrive with M2");
  const Vec<13> ref = c.exact(c.t_final);
  std::vector<ConvergenceRow> rows;
  for (int k = 0; k < levels; ++k) {
    const double h = h0 / std::ldexp(1.0, k);
    const RunResult r = integrate_case(c, integ, h);
    rows.push_back({c.name, integ, h, state_error(r.x_final, ref), r.nfev});
  }
  return rows;
}

OrderFit observed_order(const std::vector<double>& h, const std::vector<double>& err, double scale) {
  const OrderFit none{std::numeric_limits<double>::quiet_NaN(), -1, -1};
  if (h.size() != err.size() || h.size() < 2) return none;
  const double floor_ = 1e3 * std::numeric_limits<double>::epsilon() * scale;

  // Finest level still above the roundoff floor (h decreases with the index).
  int last = -1;
  for (std::size_t i = 0; i < err.size(); ++i)
    if (err[i] > floor_) last = static_cast<int>(i);
  if (last < 1) return none;

  // Local orders p_k between level k and k+1; extend the window backwards from the finest level while they
  // agree.
  auto p = [&](int k) {
    const auto i = static_cast<std::size_t>(k);
    return std::log(err[i] / err[i + 1]) / std::log(h[i] / h[i + 1]);
  };
  int first = last - 1;
  while (first > 0 && err[static_cast<std::size_t>(first - 1)] > floor_ &&
         std::abs(p(first - 1) - p(first)) <= 0.1 * std::abs(p(first)))
    --first;

  // Least-squares slope of log(err) against log(h) over [first, last].
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  const double m = last - first + 1;
  for (int k = first; k <= last; ++k) {
    const double x = std::log(h[static_cast<std::size_t>(k)]), y = std::log(err[static_cast<std::size_t>(k)]);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
  }
  return {(m * sxy - sx * sy) / (m * sxx - sx * sx), first, last};
}

}  // namespace flightsim::validation
