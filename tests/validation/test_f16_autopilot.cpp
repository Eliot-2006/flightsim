// PRD 8.7 control validation (C1-C4) from the COMMITTED gains and scenarios, plus the closed-loop energy
// budget. C5 (loop margins) is asserted in python/tests/test_design.py; here the nonlinear model is only
// observed.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/sim/aircraft_runner.hpp"

using namespace flightsim;
namespace {
const std::filesystem::path root = FLIGHTSIM_SOURCE_DIR;
constexpr double kFt = units::ft;  // PRD targets are in feet; the runner reports SI

struct Run {
  sim::AircraftScenario scenario;
  sim::AircraftRunOutput out;
  sim::ControlReport report;
};
Run run(const char* name) {
  Run r{sim::load_aircraft_scenario(root / "scenarios" / name), {}, {}};
  r.out = sim::run_aircraft(r.scenario);
  r.report = sim::control_report(r.scenario, r.out);
  // Every run: finishes, and the table-domain monitor stays clean (a response from extrapolated tables is not
  // a validation result).
  EXPECT_EQ(r.out.reason, sim::StopReason::FinalTime) << name << ": " << r.out.message;
  EXPECT_FALSE(r.out.out_of_domain) << name << ": " << r.out.first_violation;
  return r;
}
void print(const char* name, const sim::StepReport& s) {
  std::cout << std::setprecision(5) << name << ' ' << s.kind << " step " << s.metrics.step / kFt
            << " ft-or-ft/s: overshoot " << s.metrics.overshoot_pct << " %, 2% settling "
            << s.metrics.settling_2pct_s << " s, steady-state |e| (last 10 s) "
            << s.metrics.steady_state_error / kFt << ", other-loop peak |e| " << s.peak_cross_deviation / kFt
            << '\n';
}
}  // namespace

TEST(F16Autopilot, C1AltitudeStep100ft) {
  const auto r = run("f16_ap_c1_alt100.toml");
  ASSERT_EQ(r.report.steps.size(), 1u);
  const auto& s = r.report.steps[0];
  print("C1", s);
  EXPECT_NEAR(s.metrics.step / kFt, 100.0, 1e-9);
  // PRD 8.7 C1: overshoot <= 10 %, 2 % settling <= 25 s, steady-state |e| <= 1 ft, airspeed deviation <= 5
  // ft/s.
  EXPECT_LE(s.metrics.overshoot_pct, 10.0);
  EXPECT_LE(s.metrics.settling_2pct_s, 25.0);
  EXPECT_FALSE(std::isnan(s.metrics.settling_2pct_s));
  // Steady state = last 10 s of a run that spans the design's 300 s linear horizon (scenario comment: the
  // altitude PI zero leaves a ~200 s tail; at 60 s the linear model itself still has 1.09 ft).
  EXPECT_LE(s.metrics.steady_state_error, 1.0 * kFt);
  EXPECT_LE(s.peak_cross_deviation, 5.0 * kFt);
}

TEST(F16Autopilot, C2AltitudeStep1000ftAntiwindupBeatsNoAntiwindup) {
  const auto with = run("f16_ap_c2_alt1000.toml"), without = run("f16_ap_c2_alt1000_nowindup.toml");
  ASSERT_EQ(with.report.steps.size(), 1u);
  ASSERT_EQ(without.report.steps.size(), 1u);
  const auto &a = with.report.steps[0].metrics, &n = without.report.steps[0].metrics;
  print("C2 anti-windup", with.report.steps[0]);
  print("C2 no anti-windup", without.report.steps[0]);
  EXPECT_TRUE(with.scenario.autopilot->antiwindup);
  EXPECT_FALSE(without.scenario.autopilot->antiwindup);
  // PRD 8.7 C2: with anti-windup overshoot <= 15 %.
  EXPECT_LE(a.overshoot_pct, 15.0);
  // M7b review acceptance: the anti-windup run must also settle (2 % within 120 s) and end within 10 ft of
  // the command (max |e| over the last 10 s of the 305 s run).
  ASSERT_FALSE(std::isnan(a.settling_2pct_s)) << "anti-windup run never settled";
  EXPECT_LE(a.settling_2pct_s, 120.0);
  EXPECT_LE(a.steady_state_error, 10.0 * kFt);
  // Anti-windup must be no worse than none on BOTH overshoot and settling time (an unsettled run counts as
  // worse).
  EXPECT_LE(a.overshoot_pct, n.overshoot_pct);
  if (!std::isnan(n.settling_2pct_s)) {
    EXPECT_LE(a.settling_2pct_s, n.settling_2pct_s);
  }
}

TEST(F16Autopilot, C3AirspeedStep50fps) {
  const auto r = run("f16_ap_c3_speed50.toml");
  ASSERT_EQ(r.report.steps.size(), 1u);
  const auto& s = r.report.steps[0];
  print("C3", s);
  EXPECT_NEAR(s.metrics.step / kFt, 50.0, 1e-9);
  // PRD 8.7 C3: overshoot <= 10 %, altitude deviation <= 20 ft.
  EXPECT_LE(s.metrics.overshoot_pct, 10.0);
  EXPECT_LE(s.peak_cross_deviation, 20.0 * kFt);
  // No target for C3 settling, but the airspeed must settle within the run (reported above).
  EXPECT_FALSE(std::isnan(s.metrics.settling_2pct_s)) << "C3 airspeed never settled within 2 %";
}

TEST(F16Autopilot, C4GustAndHeadwindChangeRecoverWithinTwoPercentBand) {
  const auto r = run("f16_ap_c4_disturbance.toml");
  ASSERT_EQ(r.report.disturbances.size(), 2u);
  for (const auto& d : r.report.disturbances) {
    std::cout << std::setprecision(5) << "C4 " << d.kind << " [" << d.t_s << ", " << d.t_end_s
              << "] s: peak altitude deviation " << d.altitude.peak_deviation / kFt << " ft, recovery "
              << d.altitude.recovery_2pct_s << " s; peak airspeed deviation "
              << d.airspeed.peak_deviation / kFt << " ft/s, recovery " << d.airspeed.recovery_2pct_s
              << " s\n";
    // PRD 8.7 C4: peak deviation and recovery time reported; altitude "returns within a 2 % band" (of its
    // peak excursion) before the next disturbance / the end of the run.
    EXPECT_GT(d.altitude.peak_deviation, 0.0) << d.kind;
    EXPECT_FALSE(std::isnan(d.altitude.recovery_2pct_s)) << d.kind << " altitude did not return to the band";
  }
  EXPECT_EQ(r.report.disturbances[0].kind, "wind_step");
  EXPECT_EQ(r.report.disturbances[1].kind, "gust");
  // The headwind change is a 20 kt step in air-relative speed (inertial velocity is continuous).
  EXPECT_NEAR(r.report.disturbances[0].airspeed.peak_deviation, 20 * 1852.0 / 3600, 1e-3);
}

// B1 with the controller in the loop (PRD 8.4, M4 amendment: the open-loop 60 s run leaves the alpha table
// domain, so the 60 s budget check is closed-loop). Augmented state [x17 | W_nc], W_nc' = (C_NB F_B).v_N +
// M_B.w.
//
// Why the balance still holds exactly per step: the controller is OUTSIDE the derivative. Within one
// integration step the commanded inputs are constants (sample and hold, dt divides the controller period), so
// the augmented system integrated by RK4 is the same smooth ODE as in the open-loop B1 run; the controller
// only changes the constants BETWEEN steps, while E and W_nc are functions of the (continuous) state. d/dt(E
// - W_nc) = 0 holds for ANY input history, so feedback cannot create a budget residual; it only moves the
// trajectory, which stays in the table domain (asserted).
namespace {
double energy(const F16Model& m, const Vec<17>& x) {
  const auto mp = m.mass_props();
  const Vec3 w = rb::omega(x);
  return 0.5 * mp.mass * rb::vel(x).squaredNorm() + 0.5 * w.dot(mp.inertia * w) - mp.mass * m.grav.g0 * x(2);
}
}  // namespace

TEST(F16AutopilotB1, ClosedLoopWorkEnergyBalance60s) {
  auto s = sim::load_aircraft_scenario(root / "scenarios/f16_ap_c1_alt100.toml");
  ASSERT_TRUE(s.autopilot.has_value());
  s.run->t_final = 60.0;  // PRD: 60 s of C1 (the +100 ft command at 5 s and the response to it)
  F16Model m = s.model;
  LookupMonitor monitor;
  m.monitor = &monitor;
  const auto tr = trim_level(m, s.trim);
  ASSERT_TRUE(tr.ok) << tr.message;

  const double dt = s.run->dt;
  const long n = std::lround(s.run->t_final / dt);
  const auto measure = [&](double t, const Vec<17>& x) {
    const Vec<13> c = x13_from_x17(m, x, t);
    return ApMeasurement{c(11), c(0), c(4), c(7)};
  };
  Autopilot ap(s.autopilot->gains, s.autopilot->antiwindup);
  Vec<18> x = Vec<18>::Zero();
  x.head<17>() = tr.x;
  const ApMeasurement m0 = measure(0, tr.x);
  ap.engage(m0, tr.u);
  ApCommand cmd{m0.h, m0.V};
  F16Input u = tr.u;
  const double E0 = energy(m, tr.x);
  double max_delta = 0, max_residual = 0;
  for (long k = 0; k < n; ++k) {
    const double t = static_cast<double>(k) * dt;
    if (k % s.autopilot->steps_per_sample == 0) {
      for (const auto& c : s.commands) {
        if (k >= std::lround(c.t_s / dt)) {
          (c.kind == sim::Command::Altitude ? cmd.h : cmd.V) = c.value;
        }
      }
      u = ap.update(measure(t, x.head<17>()), cmd).u;
    }
    const auto f = [&](double time, const Vec<18>& state) -> Vec<18> {
      const Vec<17> a = state.head<17>();
      Vec<18> dx = Vec<18>::Zero();
      dx.head<17>() = m.deriv(time, a, u);
      const auto fm = m.evaluate(time, a, u).fm;
      dx(17) = quat_rotate(rb::quat(a), fm.F_B).dot(rb::vel(a)) + fm.M_B.dot(rb::omega(a));
      return dx;
    };
    x = rk4_step<18>(f, t, x, dt);
    renormalize_quat(x);
    ASSERT_TRUE(x.allFinite()) << "t=" << t + dt;
    const double delta = energy(m, x.head<17>()) - E0;
    max_delta = std::max(max_delta, std::abs(delta));
    max_residual = std::max(max_residual, std::abs(delta - x(17)));
  }
  EXPECT_FALSE(monitor.out_of_domain) << monitor.first_violation;
  ASSERT_GT(max_delta, 1e5) << "the +100 ft maneuver must exchange energy (about m g dh = 2.8 MJ)";
  const double ratio = max_residual / max_delta;
  std::cout << std::setprecision(6) << "B1 closed loop (C1, 60 s, dt=0.01): max|dE| = " << max_delta
            << " J, max|dE - W_nc| = " << max_residual << " J, ratio = " << ratio << '\n';
  // PRD 8.4 B1 / M4 amendment: |dE - W_nc| / max|dE| < 1e-6 at dt = 0.01.
  EXPECT_LT(ratio, 1e-6);

  // The hand-rolled loop above must be the runner's loop: same final state as run_aircraft on the same
  // scenario.
  const auto reference = sim::run_aircraft(s);
  ASSERT_EQ(reference.reason, sim::StopReason::FinalTime);
  const Vec<17> a = x.head<17>(), b = reference.x.back();
  // 1e-9 relative: identical arithmetic per component (fp-contract off); slack only for compiler differences.
  EXPECT_LT((a - b).norm(), 1e-9 * (1 + b.norm()));
}
