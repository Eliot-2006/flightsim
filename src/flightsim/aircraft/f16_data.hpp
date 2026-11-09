#pragma once

#include <array>
#include <filesystem>

#include "flightsim/math/interp.hpp"

namespace flightsim {

// Per-surface actuator parameters. TOML table: [actuators.elevator|aileron|rudder].
// Positive deflection conventions (S&L, see data/f16/SOURCES.md): elevator trailing-edge-down, aileron
// left-roll-command (right TE down), rudder left-yaw-command (TE left).
struct F16Actuator {
  double pos_min;     // rad.  TOML position_limits[0]  (deg -> rad)
  double pos_max;     // rad.  TOML position_limits[1]  (deg -> rad)
  double rate_limit;  // rad/s. TOML rate_limit         (deg/s -> rad/s)
  double lag_rate;    // 1/s.   TOML lag_rate, first-order lag inverse time constant; unchanged
};

// S&L `adc` atmosphere constants, kept in the ORIGINAL IMPERIAL units of the TOML (not converted).
// They exist so a test can pin env/atmosphere.cpp against the data file; flight code uses atmosphere().
// TOML table: [environment.*], one scalar table per field named below.
struct F16SlAdc {
  double rho0;              // slug/ft^3. environment.rho0
  double temperature0;      // degR.      environment.temperature0
  double temperature_high;  // degR.      environment.temperature_high (used at/above break_altitude)
  double break_altitude;    // ft.        environment.temperature_break_altitude
  double altitude_factor;   // 1/ft.      environment.temperature_altitude_factor (tfac = 1 - factor*h_ft)
  double density_exponent;  // 1.         environment.density_exponent
  double gamma;             // 1.         environment.gamma_air
  double gas_constant;      // ft^2/(s^2*degR). environment.gas_constant
};

// Everything the Stevens & Lewis low-fidelity F-16 model needs, SI, angles in rad.
// Loaded by load_f16_data() from data/f16/f16_sl.toml. Coefficient VALUES are never modified; only axes and
// physical constants are unit-converted. Frames: body FRD. Altitude is geometric, positive up.
//
// Table argument order / axes (call as table(first, second)); the lookup extrapolates linearly from the
// boundary cell without limit -- pass a LookupMonitor to detect that (see interp.hpp).
//
// Alpha, beta, elevator inputs to every table are in RADIANS (the loader converted the axes).
// Aero tables are dimensionless coefficients. Sign conventions of the data file apply (SOURCES.md):
// Cl and Cn tables are defined for |beta| and odd in beta -> evaluate with |beta| and multiply by sign(beta);
// the DLDA/DLDR/DNDA/DNDR effectiveness tables use SIGNED beta and must NOT be sign-folded.
struct F16Data {
  // ---- [mass_geometry.*]
  // -------------------------------------------------------------------------------------
  double weight_nominal;  // N.      NASA nominal weight W (20500 lbf). Informational only; mass below is what
                          // S&L uses.
                          //         TOML mass_geometry.W
  double mass;  // kg.     = (1 / reciprocal_mass) * slug. S&L uses RM = 0.00157 1/slug, NOT W/g (differs by
                //         ~1e-4 relative; do not "fix"). TOML mass_geometry.reciprocal_mass
  double Ixx, Iyy, Izz;  // kg*m^2. Body-axis moments of inertia. TOML mass_geometry.Ixx / Iyy / Izz
  double Ixz;            // kg*m^2. POSITIVE product of inertia int(x z dm) as printed by S&L (982 slug ft^2).
                         //         The FRD inertia tensor off-diagonal elements I(0,2)=I(2,0) are -Ixz.
                         //         TOML mass_geometry.Ixz
  double S;              // m^2.    Wing reference area. TOML mass_geometry.S
  double b;              // m.      Wing span. TOML mass_geometry.b
  double cbar;           // m.      Mean aerodynamic chord. TOML mass_geometry.cbar
  double
      xcg_ref;   // 1.      Reference CG as a FRACTION of cbar (0.35, not percent). TOML mass_geometry.xcg_ref
  double hx;     // kg*m^2/s. Engine rotor angular momentum along +x body. TOML mass_geometry.hx
  double g_sl;   // m/s^2.  Gravity used by the S&L reference (32.17 ft/s^2 = 9.80616..., not units::g0).
                 //         TOML environment.g
  double rtod;   // deg/rad, UNCONVERTED. S&L's rounded 57.29578, distinct from 180/pi; only needed to
                 //         reproduce S&L check-case numbers exactly. TOML environment.rtod
  F16SlAdc adc;  // see F16SlAdc

  // ---- [aero.*] static force/moment coefficient tables
  // ---------------------------------------------------------
  Table2D cx;    // Cx(alpha, elevator).           axes: alpha [rad], elevator [rad].  TOML aero.cx
  Table2D cm;    // Cm(alpha, elevator).           axes: alpha, elevator.              TOML aero.cm
  Table2D cl;    // Cl(alpha, |beta|)  per |beta|. axes: alpha, abs_beta [rad, >= 0]. TOML aero.cl
  Table2D cn;    // Cn(alpha, |beta|)  per |beta|. axes: alpha, abs_beta.              TOML aero.cn
  Table2D dlda;  // dCl for aileron at normalized aileron=1 (delta_a / aileron_norm).  axes: alpha, SIGNED
                 // beta.  aero.dlda
  Table2D dldr;  // dCl for rudder at normalized rudder=1 (delta_r / rudder_norm).      axes: alpha, SIGNED
                 // beta.  aero.dldr
  Table2D dnda;  // dCn for aileron, as dlda.                                           aero.dnda
  Table2D dndr;  // dCn for rudder, as dldr.                                            aero.dndr
  Table1D cz;    // Cz base(alpha). axis: alpha [rad]. S&L: Cz = cz(alpha)*(1 - (beta/cz_beta_scale)^2)
               //         + cz_elevator_gain*(elevator/cz_elevator_scale). TOML aero.cz (values + scalar keys
               //         below)
  double cz_beta_scale;      // rad. 57.3 deg -> rad  (aero.cz.beta_scale; NOT the same constant as rtod)
  double cz_elevator_gain;   // 1.   -0.19            (aero.cz.elevator_gain)
  double cz_elevator_scale;  // rad. 25 deg -> rad    (aero.cz.elevator_scale)
  // Side force: Cy = cy_beta_gain*beta + cy_aileron_gain*(aileron/aileron_norm) +
  // cy_rudder_gain*(rudder/rudder_norm)
  double cy_beta_gain;     // 1/rad. TOML -0.02 1/deg divided by units::deg  (aero.cy.beta_gain)
  double cy_aileron_gain;  // 1.     (aero.cy.aileron_gain)
  double cy_rudder_gain;   // 1.     (aero.cy.rudder_gain)
  double
      aileron_norm;  // rad.   20 deg -> rad; aileron deflection divisor (aero.control_normalization.aileron)
  double rudder_norm;  // rad.   30 deg -> rad; rudder deflection divisor  (aero.control_normalization.rudder)

  // ---- [aero.damping] dynamic derivatives, each a function of alpha [rad] (TOML columns of
  // aero.damping.values) ---- Dimensionless. Non-dimensional rates: CXq, CZq, Cmq multiply q*cbar/(2V); CYr,
  // Clr, Cnr multiply r*b/(2V); CYp, Clp, Cnp multiply p*b/(2V).
  Table1D cxq, cyr, cyp, czq, clr, clp, cmq, cnr, cnp;

  // ---- [engine.idle|military|maximum] thrust tables
  // ---------------------------------------------------------------- Thrust(altitude, Mach). Axis 0:
  // geometric altitude [m] (TOML ft -> m); axis 1: Mach, dimensionless. Output: N along +x body (TOML lbf ->
  // N). Idle can be negative (drag) -- intentional data.
  Table2D thrust_idle;  // power level 0.    TOML engine.idle
  Table2D thrust_mil;   // power level 50.   TOML engine.military
  Table2D thrust_max;   // power level 100.  TOML engine.maximum
  std::array<double, 3>
      power_levels;  // percent {0, 50, 100}: levels of idle/mil/max. engine.thrust_power_levels.values

  // ---- Engine constants, percent / 1/s, UNCHANGED from the TOML
  // --------------------------------------------------------- Throttle gearing (throttle in [0,1] ->
  // commanded power in percent): throttle <= gear_breakpoint uses the low branch:
  //   power_cmd = gear_low_slope * throttle ; otherwise power_cmd = gear_high_slope * throttle +
  //   gear_high_intercept.
  double gear_breakpoint;  // 1.       engine.throttle_gearing.breakpoint
  double gear_low_slope;   // percent  engine.throttle_gearing.low_slope
  double gear_high_slope;  // percent  engine.throttle_gearing.high_slope
  double
      gear_high_intercept;  // percent  engine.throttle_gearing.high_intercept (retained as printed, -117.38)
  // rtau = inverse engine time constant [1/s] as a function of delta_power = |power_cmd - power| [percent]:
  //   delta_power <= rtau_dp_low -> rtau_low ; >= rtau_dp_high -> rtau_high ; between -> mid_intercept +
  //   mid_slope*dp.
  double rtau_dp_low;         // percent    engine.rtau.delta_power[0]
  double rtau_dp_high;        // percent    engine.rtau.delta_power[1]
  double rtau_low;            // 1/s        engine.rtau.values[0]
  double rtau_high;           // 1/s        engine.rtau.values[1]
  double rtau_mid_intercept;  // 1/s        engine.rtau.middle_intercept
  double rtau_mid_slope;      // 1/(s*percent) engine.rtau.middle_slope
  // Power lag (see data/f16/SOURCES.md "Engine"): when power and command straddle power_threshold the command
  // is replaced by power_up_target (accelerating) / power_down_target (decelerating); at power >=
  // power_threshold the rate is power_high_rate instead of rtau.
  double power_threshold;    // percent    engine.power_lag.power_threshold
  double power_up_target;    // percent    engine.power_lag.transition_up_target
  double power_down_target;  // percent    engine.power_lag.transition_down_target
  double power_high_rate;    // 1/s        engine.power_lag.high_power_rate

  // ---- [actuators.*]
  // ------------------------------------------------------------------------------------------------
  F16Actuator elevator;  // actuators.elevator
  F16Actuator aileron;   // actuators.aileron
  F16Actuator rudder;    // actuators.rudder
  double throttle_min;   // 1. actuators.throttle.position_limits[0]
  double throttle_max;   // 1. actuators.throttle.position_limits[1]
  // [actuators.aerobench] only records that AeroBench has no actuator dynamics; it is validated present, not
  // stored.
};

// Loads and unit-converts the S&L data file. Throws std::runtime_error whose message names the file and the
// offending TOML key/table path when: the file is unreadable or malformed; schema_version != 1; a required
// table/key is missing; an unknown table (a table with a `source` key that the loader does not consume) is
// present; a units/axes label differs from what the loader converts; a table fails the Table1D/Table2D
// checks.
F16Data load_f16_data(const std::filesystem::path& toml_path);

}  // namespace flightsim
