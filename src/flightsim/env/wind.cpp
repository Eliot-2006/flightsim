#include "flightsim/env/wind.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace flightsim {

Vec3 WindModel::wind_N(double t) const {
  if (!std::isfinite(t) || !constant_N.allFinite()) {
    throw std::invalid_argument("wind: time and constant wind must be finite");
  }
  if (!gust) return constant_N;
  const auto& g = *gust;
  const double norm = g.direction_N.stableNorm();
  if (!std::isfinite(g.t_start_s) || !std::isfinite(g.duration_s) || g.duration_s <= 0 ||
      !std::isfinite(g.amplitude_mps) || !g.direction_N.allFinite() || !std::isfinite(norm) || norm == 0) {
    throw std::invalid_argument("wind: gust requires finite parameters, positive duration and direction");
  }
  const double elapsed = t - g.t_start_s;
  if (elapsed <= 0 || elapsed >= g.duration_s) return constant_N;
  return constant_N +
         (0.5 * g.amplitude_mps * (1 - std::cos(2 * std::numbers::pi * elapsed / g.duration_s))) *
             (g.direction_N / norm);
}

}  // namespace flightsim
