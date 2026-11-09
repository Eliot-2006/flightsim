#include "flightsim/analysis/linearize.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/analysis/jacobian.hpp"
#include "flightsim/core/integrators.hpp"

namespace flightsim {
LinearModel linearize(const F16Model& m, const TrimResult& tr, double eps) {
  if (!tr.ok || m.actuators) {
    throw std::invalid_argument("linearize requires a successful trim and ideal surfaces");
  }
  LinearModel lm;
  lm.x0 = x13_from_x17(m, tr.x);
  lm.u0 << tr.u.throttle, tr.u.de, tr.u.da, tr.u.dr;
  lm.A = jacobian_central<13, 13>([&](const Vec<13>& x) { return f13(m, x, tr.u); }, lm.x0, Vec<13>::Ones(),
                                  eps);
  lm.B = jacobian_central<13, 4>([&](const Vec<4>& u) { return f13(m, lm.x0, {u(0), u(1), u(2), u(3)}); },
                                 lm.u0, Vec<4>::Ones(), eps);
  return lm;
}
std::vector<SweepPoint> linearization_sweep(const F16Model& m, const TrimResult& tr) {
  std::vector<SweepPoint> points;
  auto previous = linearize(m, tr, 1e-1);
  for (int i = 2; i <= 10; ++i) {
    const double eps = std::pow(10.0, -i);
    const auto lm = linearize(m, tr, eps);
    const double change = std::max((lm.A - previous.A).cwiseAbs().maxCoeff() / lm.A.cwiseAbs().maxCoeff(),
                                   (lm.B - previous.B).cwiseAbs().maxCoeff() / lm.B.cwiseAbs().maxCoeff());
    points.push_back({eps, change});
    previous = lm;
  }
  return points;
}
ElevatorResponse elevator_response(const F16Model& model, const TrimResult& tr, const LinearModel& lm,
                                   double step) {
  F16Model m = model;
  LookupMonitor monitor;
  monitor.record_cells = true;
  m.monitor = &monitor;
  F16Input u = tr.u;
  u.de += step;
  ElevatorResponse result;
  Vec<17> x = tr.x;
  Vec<13> dx = Vec<13>::Zero();
  Vec<4> du = Vec<4>::Zero();
  du(1) = step;
  m.evaluate(0, x, u);
  const auto cells = monitor.cells;
  monitor.reset();
  auto f = [&](double t, const Vec<17>& y) -> Vec<17> {
    monitor.cells.clear();
    const Vec<17> dy = m.deriv(t, y, u);
    if (monitor.cells != cells) {
      if (!result.cells_changed) {
        result.first_cell_change_time = t;
      }
      result.cells_changed = true;
    }
    return dy;
  };
  for (int i = 0; i <= 500; ++i) {
    const Vec<13> nonlinear = x13_from_x17(m, x) - lm.x0;
    // Position integrators have a nonzero equilibrium velocity; remove its nominal trajectory.
    Vec<13> perturbation = nonlinear;
    perturbation(9) -= tr.x(rb::kVel) * i * .01;
    perturbation(10) -= tr.x(rb::kVel + 1) * i * .01;
    perturbation(11) += tr.x(rb::kVel + 2) * i * .01;
    const double error = (perturbation.head<9>() - dx.head<9>()).norm();
    result.max_error = std::max(result.max_error, error);
    result.points.push_back({i * .01, perturbation, dx});
    if (i == 500) {
      break;
    }
    x = rk4_step<17>(f, i * .01, x, .01);
    rb::quat(x) = quat_normalized(rb::quat(x));
    dx = rk4_step<13>([&](double, const Vec<13>& y) -> Vec<13> { return lm.A * y + lm.B * du; }, i * .01, dx,
                      .01);
  }
  result.out_of_domain = monitor.out_of_domain;
  result.first_violation = monitor.first_violation;
  return result;
}
}  // namespace flightsim
