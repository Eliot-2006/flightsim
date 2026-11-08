#pragma once

#include <numbers>

namespace flightsim::units {

// Imperial <-> SI conversion constants, defined once (PRD 4.3). Multiply an imperial value to get SI.
inline constexpr double ft = 0.3048;                     // m, exact
inline constexpr double slug = 14.59390294;              // kg
inline constexpr double lbf = 4.4482216152605;           // N
inline constexpr double deg = std::numbers::pi / 180.0;  // rad
inline constexpr double g0 = 9.80665;                    // m/s^2, standard gravity (exact)

}  // namespace flightsim::units
