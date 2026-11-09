#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;

namespace {

F16Model model() {
  F16Model m{load_f16_data(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml")};
  m.grav.g0 = m.data.g_sl;
  return m;
}

F16Input controls(double t, bool lateral_only = false, double sign = 1) {
  F16Input u{0.1678, -0.6531 * units::deg, 0, 0};
  if (!lateral_only) {
    if (t >= 1 && t < 2) u.de += units::deg;
    if (t >= 2 && t < 3) u.de -= units::deg;
  }
  if (t >= 5 && t < 6) u.da = sign * 2 * units::deg;
  if (t >= 6 && t < 7) u.da = -sign * 2 * units::deg;
  return u;
}

Vec<17> initial(const F16Model& m) {
  Vec<13> s = Vec<13>::Zero();
  s(0) = 152.4;
  s(1) = s(4) = 0.0594;
  s(11) = 3048;
  s(12) = 64.94 * 0.1678;
  return x17_from_x13(m, s, controls(0));
}

double energy(const F16Model& m, const Vec<17>& x) {
  const auto mp = m.mass_props();
  const Vec3 w = rb::omega(x);
  return 0.5 * mp.mass * rb::vel(x).squaredNorm() + 0.5 * w.dot(mp.inertia * w) - mp.mass * m.grav.g0 * x(2);
}

Vec3 momentum(const F16Model& m, const Vec<17>& x) {
  const auto mp = m.mass_props();
  return quat_rotate(rb::quat(x), mp.inertia * rb::omega(x) + mp.h_rotor);
}

struct Run {
  std::vector<Vec<17>> states;
  double residual = 0;
  bool out_of_domain = false;
  std::string first_violation;
};

// N=18 adds W_nc; N=20 adds J_N. Inputs are held constant across all four stages per PRD 6.2.
template <int N = 17>
Run run(F16Model m, double dt, double duration, bool lateral_only = false, double sign = 1) {
  LookupMonitor monitor;
  m.monitor = &monitor;
  Vec<N> x = Vec<N>::Zero();
  x.template head<17>() = initial(m);
  const double E0 = energy(m, x.template head<17>());
  const Vec3 H0 = momentum(m, x.template head<17>());
  double max_delta = 0, max_residual = 0;
  const int steps = static_cast<int>(std::lround(duration / dt));
  Run result;
  result.states.reserve(static_cast<std::size_t>(steps + 1));
  result.states.push_back(x.template head<17>());
  for (int i = 0; i < steps; ++i) {
    const double t = static_cast<double>(i) * dt;
    // Every event is aligned to all dt values in this file; midpoint avoids binary endpoint ambiguity.
    const F16Input u = controls(t + 0.5 * dt, lateral_only, sign);
    auto f = [&](double time, const Vec<N>& state) -> Vec<N> {
      const Vec<17> aircraft = state.template head<17>();
      Vec<N> dx = Vec<N>::Zero();
      const bool was_outside = monitor.out_of_domain;
      dx.template head<17>() = m.deriv(time, aircraft, u);
      if (!was_outside && monitor.out_of_domain) {
        std::cout << std::setprecision(12) << "First out-of-domain lookup dt=" << dt << " t=" << time << " "
                  << monitor.first_violation << '\n';
      }
      if constexpr (N == 18 || N == 20) {
        const auto fm = m.evaluate(time, aircraft, u).fm;
        if constexpr (N == 18) {
          dx(17) = quat_rotate(rb::quat(aircraft), fm.F_B).dot(rb::vel(aircraft)) +
                   fm.M_B.dot(rb::omega(aircraft));
        } else {
          dx.template tail<3>() = quat_rotate(rb::quat(aircraft), fm.M_B);
        }
      }
      return dx;
    };
    x = rk4_step<N>(f, t, x, dt);
    renormalize_quat(x);
    if (!x.allFinite()) {
      ADD_FAILURE() << "non-finite aircraft state at t=" << t + dt;
      break;
    }
    result.states.push_back(x.template head<17>());
    if constexpr (N == 18) {
      const double delta = energy(m, x.template head<17>()) - E0;
      max_delta = std::max(max_delta, std::abs(delta));
      max_residual = std::max(max_residual, std::abs(delta - x(17)));
    } else if constexpr (N == 20) {
      const Vec3 delta = momentum(m, x.template head<17>()) - H0;
      max_delta = std::max(max_delta, delta.norm());
      max_residual = std::max(max_residual, (delta - x.template tail<3>()).norm());
    }
  }
  result.residual = max_delta > 0 ? max_residual / max_delta : 0;
  result.out_of_domain = monitor.out_of_domain;
  result.first_violation = monitor.first_violation;
  EXPECT_FALSE(monitor.out_of_domain) << monitor.first_violation;
  EXPECT_EQ(result.states.size(), static_cast<std::size_t>(steps + 1));
  return result;
}

Table2D beta_even(const Table2D& table) {
  Eigen::MatrixXd values(static_cast<Eigen::Index>(table.x().size()),
                         static_cast<Eigen::Index>(table.y().size()));
  for (std::size_t i = 0; i < table.x().size(); ++i) {
    for (std::size_t j = 0; j < table.y().size(); ++j) {
      values(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          0.5 * (table(table.x()[i], table.y()[j]) + table(table.x()[i], -table.y()[j]));
    }
  }
  return Table2D(table.name(), table.x(), table.y(), values);
}

Vec<13> mirror_asymmetry(const F16Model& m) {
  constexpr double dt = 0.01;
  const auto positive = run(m, dt, 10, true, 1);
  const auto negative = run(m, dt, 10, true, -1);
  Vec<13> max_difference = Vec<13>::Zero(), scale = Vec<13>::Zero();
  for (std::size_t i = 0; i < positive.states.size(); ++i) {
    const double t = static_cast<double>(i) * dt;
    const auto a = x13_from_x17(m, positive.states[i], t);
    Vec<13> b = x13_from_x17(m, negative.states[i], t);
    // beta, phi, psi, p, r, East are odd; the other chart states are even.
    for (int j : {2, 3, 5, 6, 8, 10}) b(j) = -b(j);
    max_difference = max_difference.cwiseMax((a - b).cwiseAbs());
    scale = scale.cwiseMax(a.cwiseAbs()).cwiseMax(b.cwiseAbs());
  }
  // Relative to each state's trajectory maximum with a scale-one floor for zero/near-zero states.
  return max_difference.cwiseQuotient(scale.cwiseMax(1.0));
}

}  // namespace

TEST(F16Consistency, B1WorkEnergy) {
  std::array<double, 3> residuals;
  const std::array<double, 3> steps = {0.04, 0.02, 0.01};
  for (std::size_t i = 0; i < steps.size(); ++i) {
    residuals[i] = run<18>(model(), steps[i], 20).residual;
    std::cout << std::setprecision(12) << "B1 dt=" << steps[i] << " residual=" << residuals[i] << '\n';
  }
  // Amended PRD 8.4 / Task 4.6: 20 s in-domain balance, relative residual <1e-6 at dt=0.01.
  EXPECT_LT(residuals[2], 1e-6);
  EXPECT_LT(residuals[1], residuals[0]);
  EXPECT_LT(residuals[2], residuals[1]);
  std::cout << "B1 observed_orders=" << std::log2(residuals[0] / residuals[1]) << ','
            << std::log2(residuals[1] / residuals[2]) << '\n';
}

TEST(F16Consistency, B2AngularImpulse) {
  const auto result = run<20>(model(), 0.01, 20);
  std::cout << std::setprecision(12) << "B2 dt=0.01 residual=" << result.residual << '\n';
  // PRD 8.4 / Task 4.6: full inertial-vector balance including rotor, <1e-6 relative.
  EXPECT_LT(result.residual, 1e-6);
}

TEST(F16Consistency, S1MirrorSymmetrized) {
  auto m = model();
  m.include_rotor = false;
  m.data.dlda = beta_even(m.data.dlda);
  m.data.dldr = beta_even(m.data.dldr);
  m.data.dnda = beta_even(m.data.dnda);
  m.data.dndr = beta_even(m.data.dndr);
  const auto asymmetry = mirror_asymmetry(m);
  std::cout << std::setprecision(12) << "S1 symmetrized_relative=" << asymmetry.transpose() << '\n';
  // PRD 8.4 / Task 4.6: mirror parity to 1e-12 relative, using each state's maximum over 10 s.
  EXPECT_LE(asymmetry.maxCoeff(), 1e-12);
}

TEST(F16Consistency, S1ProductionAsymmetryReported) {
  const auto asymmetry = mirror_asymmetry(model());
  std::cout << std::setprecision(12) << "S1 production_relative=" << asymmetry.transpose()
            << " max=" << asymmetry.maxCoeff() << '\n';
  // Deliberately no symmetry assertion: published effectiveness tables and the rotor break parity.
}

TEST(F16Consistency, W1GalileanInvariance) {
  const auto still = model();
  auto windy = still;
  windy.wind.constant_N = Vec3(10, -5, 0);
  constexpr double dt = 0.01;
  const auto a = run(still, dt, 20);
  const auto b = run(windy, dt, 20);
  Vec<13> difference = Vec<13>::Zero(), scale = Vec<13>::Zero();
  double max_track_error = 0;
  for (std::size_t i = 0; i < a.states.size(); ++i) {
    const double t = static_cast<double>(i) * dt;
    const auto xa = x13_from_x17(still, a.states[i], t);
    const auto xb = x13_from_x17(windy, b.states[i], t);
    difference = difference.cwiseMax((xa - xb).cwiseAbs());
    scale = scale.cwiseMax(xa.cwiseAbs()).cwiseMax(xb.cwiseAbs());
    const Vec3 track_error = rb::pos(b.states[i]) - rb::pos(a.states[i]) - windy.wind.constant_N * t;
    max_track_error = std::max(max_track_error, track_error.norm());
  }
  // PRD 8.4 / Task 4.6: air-relative states <=1e-9 relative; position offset <=1e-6 m.
  for (int i : {0, 1, 2, 3, 4, 5, 6, 7, 8, 11}) {
    EXPECT_LE(difference(i) / std::max(1.0, scale(i)), 1e-9) << i;
  }
  EXPECT_LE(max_track_error, 1e-6);
  std::cout << std::setprecision(12)
            << "W1 relative=" << difference.cwiseQuotient(scale.cwiseMax(1.0)).transpose()
            << " track_error_m=" << max_track_error << '\n';
}

TEST(F16Consistency, AircraftConvergenceMeasured) {
  const auto m = model();
  const std::vector<double> steps = {0.04, 0.02, 0.01, 0.005, 0.0025};
  // Amended PRD 6.3, non-smooth aircraft only: h_min/64 reference, h_min/32 self-convergence estimate.
  const double h_ref = steps.back() / 64;
  const auto coarse_reference = run(m, 2 * h_ref, 20).states.back();
  const auto reference = run(m, h_ref, 20).states.back();
  std::vector<Vec<17>> states;
  for (double h : steps) states.push_back(run(m, h, 20).states.back());
  const double pos_reference_error = (rb::pos(coarse_reference) - rb::pos(reference)).norm();
  const double att_reference_error = attitude_error(rb::quat(coarse_reference), rb::quat(reference));
  std::vector<double> position, attitude;
  for (std::size_t i = 0; i < steps.size(); ++i) {
    position.push_back((rb::pos(states[i]) - rb::pos(reference)).norm());
    attitude.push_back(attitude_error(rb::quat(states[i]), rb::quat(reference)));
    std::cout << std::setprecision(12) << "Aircraft dt=" << steps[i]
              << " position_error_m=" << position.back() << " attitude_error_rad=" << attitude.back() << '\n';
  }
  auto measure = [&](const char* metric, const std::vector<double>& errors, double reference_error) {
    std::vector<std::size_t> levels;
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < steps.size(); ++i) {
      // Amended PRD 6.3, non-smooth case: least-squares fit over reference-resolved levels.
      if (errors[i] > 0 && errors[i] >= 10 * reference_error) {
        const double log_h = std::log(steps[i]), log_error = std::log(errors[i]);
        sx += log_h;
        sy += log_error;
        sxx += log_h * log_h;
        sxy += log_h * log_error;
        levels.push_back(i);
      }
    }
    ASSERT_GE(levels.size(), 2U) << metric << ": insufficient reference-resolved levels";
    const double n = static_cast<double>(levels.size());
    const double order = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    std::cout << "Aircraft " << metric << "_order=" << order
              << " reference_error_estimate=" << reference_error << " reference_dt=" << h_ref
              << " levels_used=";
    for (const auto level : levels) {
      std::cout << level << "(dt=" << steps[level] << ") ";
    }
    std::cout << '\n';
    // Task 4.6 / amended PRD 6.3: measure reduced order on non-smooth tables, bounds unchanged.
    EXPECT_GE(order, 0.9) << metric;
    EXPECT_LE(order, 4.2) << metric;
  };
  measure("position", position, pos_reference_error);
  measure("attitude", attitude, att_reference_error);
}
