#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>

#include "flightsim/aircraft/f16_data.hpp"
#include "flightsim/env/atmosphere.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;

namespace {

// Rankine -> kelvin, used only to express the S&L (imperial) hand-evaluated references in SI.
constexpr double kK_per_R = 5.0 / 9.0;

struct Imperial {
  double h_ft, rho_slug_ft3, T_R, a_fts, p_lbf_ft2;
};

// SI view of an imperial S&L adc result.
AtmState to_si(const Imperial& r) {
  return {r.rho_slug_ft3 * units::slug / (units::ft * units::ft * units::ft), r.T_R * kK_per_R,
          r.p_lbf_ft2 * units::lbf / (units::ft * units::ft), r.a_fts * units::ft};
}

}  // namespace

// Reference: U.S. Standard Atmosphere, 1976 (NOAA-S/T 76-1562; NASA-TM-X-74335, NTRS 19770009539).
//   Temperature, T, p, rho:  Part 4 "Main tables", Table I "Geometric Altitude, Metric Units" (printed pp.
//   51-61). Speed of sound:          Table III "Geometric Altitude, Metric Units" (printed pp. 101-109).
// Rows are read at the stated geometric altitude Z (m). Table I prints P in millibar (x100 -> Pa) with 5
// significant digits, and Table III prints speed of sound to 0.01 m/s. Tolerance 1e-4 relative: the tables
// are rounded to 5 significant digits (<= 4.5e-5), and the Z-rows were evidently computed at geopotential H
// rounded to 1 m (Z=15000 m lists H=14965, true 14964.69), which shifts p by up to g/(R T) * 0.5 m ~ 8e-5 in
// the isothermal layer (the printed 12111 vs 12111.8 Pa at 15 km).
TEST(Atmosphere, Us76MatchesStandardTable) {
  struct Row {
    double z_m, T_K, p_mb, rho, a;
  };
  const Row rows[] = {
      {0, 288.150, 1013.25, 1.2250, 340.29},       // Table I p.52 / Table III p.101
      {1000, 281.651, 898.76, 1.1117, 336.43},     // 8.9876e2 mb, 1.1117 kg/m^3
      {5000, 255.676, 540.48, 0.73643, 320.55},    // 5.4048e2 mb, 7.3643e-1 kg/m^3 (Table I p.55 / III p.103)
      {11000, 216.774, 226.99, 0.36480, 295.15},   // 2.2699e2 mb, 3.6480e-1 (Table I p.59 / III p.107)
      {15000, 216.650, 121.11, 0.19476, 295.07},   // 1.2111e2 mb, 1.9476e-1
      {20000, 216.650, 55.293, 0.088910, 295.07},  // 5.5293e1 mb, 8.8910e-2 (Table I p.61 / III p.109)
  };
  constexpr double kTol = 1e-4;
  for (const Row& r : rows) {
    const AtmState s = atmosphere(AtmosphereKind::Us76, r.z_m);
    EXPECT_NEAR(s.T, r.T_K, kTol * r.T_K) << "Z=" << r.z_m;
    EXPECT_NEAR(s.p, r.p_mb * 100.0, kTol * r.p_mb * 100.0) << "Z=" << r.z_m;
    EXPECT_NEAR(s.rho, r.rho, kTol * r.rho) << "Z=" << r.z_m;
    EXPECT_NEAR(s.a, r.a, kTol * r.a) << "Z=" << r.z_m;
  }
}

// Hand-evaluated S&L adc: tfac = 1 - 0.703e-5 h_ft; T = 519 tfac degR (390 degR for h >= 35,000 ft);
// rho = 2.377e-3 tfac^4.14 slug/ft^3; a = sqrt(1.4 * 1716.3 * T) ft/s; p = rho * 1716.3 * T lbf/ft^2.
// Tolerance 1e-12 relative: literals carry 14+ digits; the slack covers pow/sqrt rounding and unit
// conversion.
TEST(Atmosphere, SlAdcMatchesRoutine) {
  const Imperial refs[] = {
      {0.0, 0.002377, 519.0, 1116.7200096711797, 2117.3358069},
      {10000.0, 0.0017577961215513001, 482.5143, 1076.7520653920287, 1455.7000374978372},
      {40000.0, 0.0006058799557951511, 390.0, 968.0391521007815, 405.54998957117493},
  };
  for (const Imperial& r : refs) {
    const AtmState want = to_si(r);
    const AtmState got = atmosphere(AtmosphereKind::SlAdc, r.h_ft * units::ft);
    EXPECT_NEAR(got.rho, want.rho, 1e-12 * want.rho) << r.h_ft;
    EXPECT_NEAR(got.T, want.T, 1e-12 * want.T) << r.h_ft;
    EXPECT_NEAR(got.p, want.p, 1e-12 * want.p) << r.h_ft;
    EXPECT_NEAR(got.a, want.a, 1e-12 * want.a) << r.h_ft;
  }
}

// atmosphere.cpp hard-codes the S&L constants; every constant in data/f16/f16_sl.toml [environment] must
// change the output if it differs, so evaluating the routine from the loaded file constants and comparing
// proves they match. Altitudes straddle the 35,000 ft temperature break so that the high-temperature and
// break-altitude constants are exercised. Tolerance 1e-13 relative: identical operations, differing only in
// constant folding.
TEST(Atmosphere, SlAdcConstantsMatchDataFile) {
  const auto d = load_f16_data(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml");
  const auto& c = d.adc;
  for (const double h_ft : {0.0, 5000.0, 34999.0, 35001.0, 50000.0}) {
    const double tfac = 1.0 - c.altitude_factor * h_ft;
    const double T_R = h_ft >= c.break_altitude ? c.temperature_high : c.temperature0 * tfac;
    const double rho = c.rho0 * std::pow(tfac, c.density_exponent);
    const Imperial r{h_ft, rho, T_R, std::sqrt(c.gamma * c.gas_constant * T_R), rho * c.gas_constant * T_R};
    const AtmState want = to_si(r);
    const double h_m = h_ft * units::ft;
    if (h_m > 20000.0) continue;  // atmosphere() range ends at 20 km; 50,000 ft is above it
    const AtmState got = atmosphere(AtmosphereKind::SlAdc, h_m);
    EXPECT_NEAR(got.rho, want.rho, 1e-13 * want.rho) << h_ft;
    EXPECT_NEAR(got.T, want.T, 1e-13 * want.T) << h_ft;
    EXPECT_NEAR(got.p, want.p, 1e-13 * want.p) << h_ft;
    EXPECT_NEAR(got.a, want.a, 1e-13 * want.a) << h_ft;
  }
}

// g = g0 (R/(R+h))^2 with R = 6,371,000 m, h = 10 km: hand value 9.80665 * 0.9968681509990346 (=
// (6371/6381)^2). Tolerance 1e-12 relative: double arithmetic on a 3-term product. The second check pins the
// invariant g (R+h)^2 = g0 R^2.
TEST(Gravity, GravityInverseSquare) {
  GravityModel m;
  m.R = 6371000.0;
  EXPECT_NEAR(m.gravity(10000.0), 9.775937052994681, 1e-12 * 9.775937052994681);
  EXPECT_NEAR(m.gravity(0.0), units::g0, 1e-15);  // 1e-15 abs: exactly g0 up to one rounding
  const double r = m.R + 25000.0;
  EXPECT_NEAR(m.gravity(25000.0) * r * r, units::g0 * m.R * m.R, 1e-12 * units::g0 * m.R * m.R);
}

TEST(Gravity, GravityConstantWhenRInfinite) {
  GravityModel m;  // defaults: g0 = units::g0, R = infinity
  EXPECT_TRUE(std::isinf(m.R));
  EXPECT_EQ(m.gravity(0.0), units::g0);
  EXPECT_EQ(m.gravity(10000.0), units::g0);
  EXPECT_EQ(m.gravity(-300.0), units::g0);
  m.g0 = 9.0;  // the field is honoured
  EXPECT_EQ(m.gravity(5000.0), 9.0);
}

TEST(Atmosphere, ThrowsOutsideRange) {
  for (const auto kind : {AtmosphereKind::SlAdc, AtmosphereKind::Us76}) {
    EXPECT_THROW(atmosphere(kind, -500.1), std::domain_error);
    EXPECT_THROW(atmosphere(kind, 20000.1), std::domain_error);
    EXPECT_THROW(atmosphere(kind, std::nan("")), std::domain_error);
    EXPECT_NO_THROW(atmosphere(kind, -500.0));  // limits are inclusive
    EXPECT_NO_THROW(atmosphere(kind, 20000.0));
  }
}
