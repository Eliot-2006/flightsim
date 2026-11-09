#include "flightsim/analysis/linearize.hpp"

#include <algorithm>
#include <cmath>
#include <map>
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
  auto cells = monitor.cells;
  const auto& d = m.data;
  std::map<int, std::string> names;
  for (const auto* table : {&d.cx, &d.cm, &d.cl, &d.cn, &d.dlda, &d.dldr, &d.dnda, &d.dndr, &d.thrust_idle,
                            &d.thrust_mil, &d.thrust_max}) {
    names.emplace(table->id(), table->name());
  }
  for (const auto* table : {&d.cz, &d.cxq, &d.cyr, &d.cyp, &d.czq, &d.clr, &d.clp, &d.cmq, &d.cnr, &d.cnp}) {
    names.emplace(table->id(), table->name());
  }
  const auto zero_intercepts = [](const Table2D& table) {
    if (table.y().front() != 0) {
      return false;
    }
    return std::all_of(table.x().begin(), table.x().end(),
                       [&](double alpha) { return table(alpha, 0) == 0; });
  };
  const bool smooth_cl = zero_intercepts(d.cl), smooth_cn = zero_intercepts(d.cn);
  double beta = x13_from_x17(m, x)(2);
  auto record = [&](const std::string& name, int from, int to, double t, bool ignored, const char* why) {
    result.cell_changes.push_back({name, from, to, t, ignored, why});
    if (!ignored) {
      if (!result.cells_changed) {
        result.first_cell_change_time = t;
      }
      result.cells_changed = true;
    }
  };
  monitor.reset();
  auto f = [&](double t, const Vec<17>& y) -> Vec<17> {
    monitor.cells.clear();
    const Vec<17> dy = m.deriv(t, y, u);
    const double da = m.actuators ? y(ac::kDa) : u.da;
    const double dr = m.actuators ? y(ac::kDr) : u.dr;
    result.max_abs_aileron = std::max(result.max_abs_aileron, std::abs(da));
    result.max_abs_rudder = std::max(result.max_abs_rudder, std::abs(dr));
    // With ideal surfaces u is constant, so exact zeros hold throughout the entire run.
    const bool inactive = !m.actuators && da == 0 && dr == 0;
    if (monitor.cells.size() != cells.size()) {
      record("lookup sequence", -1, -1, t, false, "lookup count changed");
    } else {
      for (std::size_t i = 0; i < cells.size(); ++i) {
        const auto& before = cells[i];
        const auto& after = monitor.cells[i];
        if (before == after) {
          continue;
        }
        const bool effectiveness = after[0] == d.dlda.id() || after[0] == d.dldr.id() ||
                                   after[0] == d.dnda.id() || after[0] == d.dndr.id();
        const bool ignored = before[0] == after[0] && effectiveness && inactive;
        record(names.at(after[0]), before[1], after[1], t, ignored,
               ignored ? "da=dr=0 exactly; effectiveness multiplied by zero" : "active table cell changed");
      }
    }
    const double next_beta = x13_from_x17(m, y)(2);
    const auto sign = [](double v) { return (v > 0) - (v < 0); };
    if (sign(beta) != sign(next_beta)) {
      for (const auto* table : {&d.cl, &d.cn}) {
        const auto id = table->id();
        const auto match = [&](const auto& cell) { return cell[0] == id; };
        const auto before = std::find_if(cells.begin(), cells.end(), match);
        const auto after = std::find_if(monitor.cells.begin(), monitor.cells.end(), match);
        if (before != cells.end() && after != monitor.cells.end()) {
          const int columns = static_cast<int>(table->y().size() - 1);
          const bool smooth = (id == d.cl.id() ? smooth_cl : smooth_cn) && (*before)[1] == (*after)[1] &&
                              (*after)[1] % columns == 0;
          record(table->name(), (*before)[1], (*after)[1], t, smooth,
                 smooth ? "beta=0 fold; zero intercept at every alpha row, same first |beta| cell"
                        : "beta fold is not smooth within the same first |beta| cell");
        }
      }
    }
    beta = next_beta;
    cells = monitor.cells;
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
