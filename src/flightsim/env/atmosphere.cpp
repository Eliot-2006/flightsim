#include "flightsim/env/atmosphere.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace flightsim {

namespace {

constexpr double kMinAlt = -500.0;      // m
constexpr double kMaxAlt = 20000.0;     // m
constexpr double kK_per_R = 5.0 / 9.0;  // K per degR (exact); not in units.hpp, which this task may not edit

// Stevens & Lewis `adc`, imperial. Must equal data/f16/f16_sl.toml [environment] (pinned by
// SlAdcConstantsMatchDataFile).
constexpr double kSlRho0 = 2.377e-3;       // slug/ft^3
constexpr double kSlT0 = 519.0;            // degR
constexpr double kSlThigh = 390.0;         // degR, used at h >= kSlBreakFt
constexpr double kSlBreakFt = 35000.0;     // ft
constexpr double kSlTfacPerFt = 0.703e-5;  // 1/ft
constexpr double kSlDensityExp = 4.14;
constexpr double kSlGamma = 1.4;
constexpr double kSlGasConst = 1716.3;  // ft^2/(s^2 degR)
// S&L `adc` computes static pressure as PS = 1715 * RHO * T, a different gas constant from the 1716.3 in
// its speed of sound; both are kept verbatim. Pressure is output-only: no force, moment or derivative reads
// it.
constexpr double kSlPressureGasConst = 1715.0;  // ft^2/(s^2 degR)

AtmState sl_adc(double h_m) {
  const double h_ft = h_m / units::ft;
  const double tfac = 1.0 - kSlTfacPerFt * h_ft;
  const double T_R = h_ft >= kSlBreakFt ? kSlThigh : kSlT0 * tfac;
  const double rho =
      kSlRho0 * std::pow(tfac, kSlDensityExp);  // slug/ft^3, density keeps using tfac above the break
  return {rho * units::slug / (units::ft * units::ft * units::ft), T_R * kK_per_R,
          rho * kSlPressureGasConst * T_R * units::lbf / (units::ft * units::ft),  // S&L PS, lbf/ft^2 -> Pa
          std::sqrt(kSlGamma * kSlGasConst * T_R) * units::ft};
}

// U.S. Standard Atmosphere 1976, 0-20 km.
constexpr double kR0 = 6356766.0;    // m, geopotential reference radius
constexpr double kT0 = 288.15;       // K
constexpr double kP0 = 101325.0;     // Pa
constexpr double kRgas = 287.05287;  // J/(kg K)
constexpr double kGamma = 1.4;
constexpr double kLapse = -0.0065;  // K/m, 0-11 km
constexpr double kH11 = 11000.0;    // m, geopotential tropopause

AtmState us76(double h_m) {
  const double H = kR0 * h_m / (kR0 + h_m);
  double T, p;
  if (H <= kH11) {
    T = kT0 + kLapse * H;
    p = kP0 * std::pow(T / kT0, -units::g0 / (kRgas * kLapse));
  } else {
    T = kT0 + kLapse * kH11;
    static const double p11 = kP0 * std::pow(T / kT0, -units::g0 / (kRgas * kLapse));
    p = p11 * std::exp(-units::g0 * (H - kH11) / (kRgas * T));
  }
  return {p / (kRgas * T), T, p, std::sqrt(kGamma * kRgas * T)};
}

}  // namespace

AtmState atmosphere(AtmosphereKind kind, double h_m) {
  if (!(h_m >= kMinAlt && h_m <= kMaxAlt))  // also rejects NaN
    throw std::domain_error("atmosphere: altitude " + std::to_string(h_m) + " m outside [-500, 20000] m");
  return kind == AtmosphereKind::SlAdc ? sl_adc(h_m) : us76(h_m);
}

double GravityModel::gravity(double h_m) const {
  if (std::isinf(R)) return g0;
  const double k = R / (R + h_m);
  return g0 * k * k;
}

}  // namespace flightsim
