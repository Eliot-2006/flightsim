#pragma once

#include <functional>
#include <limits>
#include <string>

#include "flightsim/aircraft/f16_model.hpp"

namespace flightsim {
struct TrimSpec {
  double V, h, gamma, xcg;
};
struct TrimResult {
  bool ok = false;
  double alpha = 0, throttle = 0, de = 0;
  Vec<17> x = Vec<17>::Zero();
  F16Input u;
  int iterations = 0;
  double residual = std::numeric_limits<double>::infinity();
  std::string message;
};
struct TrimSolution {
  Vec3 z;
  bool ok;
  int iterations;
  double residual;
};
// Fixed three-variable LM system used by trim; domain_error rejects an infeasible trial.
TrimSolution solve_trim_system(const std::function<Vec3(const Vec3&)>& residual, Vec3 initial);
TrimResult trim_level(const F16Model&, const TrimSpec&);
}  // namespace flightsim
