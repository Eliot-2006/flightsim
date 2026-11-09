#pragma once

#include "flightsim/aircraft/f16_data.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/env/atmosphere.hpp"
#include "flightsim/env/wind.hpp"

namespace flightsim {

namespace ac {
inline constexpr int kPower = 13, kDe = 14, kDa = 15, kDr = 16, kN = 17;
}

struct F16Input {
  double throttle = 0, de = 0, da = 0, dr = 0;  // throttle [0,1]; surfaces rad, S&L signs
};

struct F16Model {
  F16Data data;
  double xcg = 0.35;
  AtmosphereKind atm = AtmosphereKind::SlAdc;
  GravityModel grav{};  // caller selects data.g_sl for S&L reproduction
  WindModel wind{};
  bool actuators = false;
  bool include_rotor = true;
  LookupMonitor* monitor = nullptr;  // diagnostics only; not shared between concurrent runs

  struct Eval {
    ForceMoment fm;
    double V, alpha, beta, mach, qbar, thrust, power_cmd;
    Vec<4> surfaces;  // effective [throttle, elevator, aileron, rudder]; surfaces in rad
  };

  // Non-gravitational aero + thrust. At V < 1e-9 m/s: alpha=beta=0, damping rates=0;
  // qbar=0 at rest, leaving engine thrust only. Out-of-table values extrapolate and flag monitor.
  // Non-finite states/inputs throw invalid_argument; atmosphere enforces its altitude domain.
  Eval evaluate(double t, const Vec<17>& x, const F16Input& u) const;
  Vec<17> deriv(double t, const Vec<17>& x, const F16Input& u) const;
  MassProps mass_props() const;  // FRD inertia off-diagonals -Ixz, rotor along +x_B
};

}  // namespace flightsim
