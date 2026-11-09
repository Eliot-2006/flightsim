#pragma once

#include <optional>

#include "flightsim/math/quat.hpp"

namespace flightsim {

struct Gust {
  double t_start_s;
  double duration_s;
  double amplitude_mps;
  Vec3 direction_N;
};

// MIL-F-8785C 1-cos discrete-gust profile, full rise-and-fall over duration_s:
// A/2 * (1-cos(2*pi*(t-t_start)/duration)), zero outside, along normalized direction_N.
// Frozen in time; spatial gradients/rotational gusts neglected (PRD 5.4).
// Invalid non-finite parameters, nonpositive duration or zero direction throw invalid_argument.
struct WindModel {
  Vec3 constant_N = Vec3::Zero();
  std::optional<Gust> gust;
  Vec3 wind_N(double t) const;
};

}  // namespace flightsim
