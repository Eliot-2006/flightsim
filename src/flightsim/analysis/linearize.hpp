#pragma once
#include <vector>

#include "flightsim/analysis/trim.hpp"

namespace flightsim {
struct LinearModel {
  Eigen::Matrix<double, 13, 13> A;
  Eigen::Matrix<double, 13, 4> B;
  Vec<13> x0;
  Vec<4> u0;
};
LinearModel linearize(const F16Model&, const TrimResult&, double eps = 6e-6);
struct SweepPoint {
  double eps, change;
};
std::vector<SweepPoint> linearization_sweep(const F16Model&, const TrimResult&);
struct ResponsePoint {
  double t;
  Vec<13> nonlinear, linear;
};
struct CellChange {
  std::string table;
  int from, to;
  double t;
  bool ignored;
  std::string reason;
};
struct ElevatorResponse {
  std::vector<ResponsePoint> points;
  double max_error = 0;
  double max_abs_aileron = 0, max_abs_rudder = 0;     // effective inputs, checked at every RK4 stage
  std::vector<CellChange> cell_changes;               // includes the two explicitly justified LN2 exclusions
  bool cells_changed = false, out_of_domain = false;  // cells_changed means an unexcluded change
  double first_cell_change_time = -1;
  std::string first_violation;
};
// 5-second RK4 response, dt=.01 s; state-scaled error norm (V in 1 m/s, others in SI).
ElevatorResponse elevator_response(const F16Model&, const TrimResult&, const LinearModel&, double de_step);
}  // namespace flightsim
