#pragma once

#include <limits>

#include "flightsim/math/units.hpp"

namespace flightsim {

enum class AtmosphereKind { SlAdc, Us76 };

struct AtmState {
  double rho, T, p, a;  // kg/m^3, K, Pa, m/s
};

// Atmosphere at GEOMETRIC altitude h_m (m, positive up). Throws std::domain_error outside [-500 m, 20 km]
// (and for NaN).
//  SlAdc: Stevens & Lewis `adc` routine -- tfac = 1 - 0.703e-5 h_ft; T = 519 tfac degR (390 degR at h >=
//  35,000 ft);
//         rho = 2.377e-3 tfac^4.14 slug/ft^3; a = sqrt(1.4 * 1716.3 * T) ft/s; p = rho * 1716.3 * T.
//         Converted to SI.
//  Us76:  U.S. Standard Atmosphere 1976, layers 0-11 km (-6.5 K/km) and 11-20 km (isothermal) in geopotential
//         altitude H = r0 h / (r0 + h), r0 = 6,356,766 m.
AtmState atmosphere(AtmosphereKind kind, double h_m);

struct GravityModel {
  double g0 = units::g0;  // m/s^2
  double R =
      std::numeric_limits<double>::infinity();  // m; body radius for the inverse-square law, inf = constant g
  double gravity(double h_m) const;             // g0 (R / (R + h))^2; R = inf -> g0
};

}  // namespace flightsim
