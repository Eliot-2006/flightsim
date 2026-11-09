#include "flightsim/analysis/trim.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/analysis/jacobian.hpp"
#include "flightsim/math/units.hpp"

namespace flightsim {
TrimSolution solve_trim_system(const std::function<Vec3(const Vec3&)>& f, Vec3 z) {
  double lambda = 1e-3;
  Vec3 r = f(z);
  for (int i = 0; i < 50; ++i) {
    if (!r.allFinite()) {
      return {z, false, i, std::numeric_limits<double>::infinity()};
    }
    // Task 5.1 fixes the scaled-residual infinity-norm threshold at 1e-10.
    if (r.lpNorm<Eigen::Infinity>() < 1e-10) {
      return {z, true, i, r.lpNorm<Eigen::Infinity>()};
    }
    const Mat3 J = jacobian_central<3, 3>(f, z, Vec3::Ones());
    const Vec3 step = (J.transpose() * J + lambda * Mat3::Identity()).ldlt().solve(-J.transpose() * r);
    Vec3 trial;
    try {
      trial = f(z + step);
    } catch (const std::domain_error&) {
      trial.setConstant(std::numeric_limits<double>::infinity());
    }
    if (trial.allFinite() && trial.squaredNorm() < r.squaredNorm()) {
      z += step;
      r = trial;
      lambda /= 10;
    } else {
      lambda *= 10;
    }
  }
  return {z, r.lpNorm<Eigen::Infinity>() < 1e-10, 50, r.lpNorm<Eigen::Infinity>()};
}

TrimResult trim_level(const F16Model& model, const TrimSpec& s) {
  TrimResult tr;
  // Trim on ideal surfaces; an actuator-enabled model then starts with its actuator states at the trimmed
  // positions (x17_from_x13 sets them from u), which makes their derivatives exactly zero.
  F16Model ideal = model;
  ideal.actuators = false;
  tr.message = "trim outside envelope or failed to converge";
  if (!std::isfinite(s.V) || s.V <= 0 || !std::isfinite(s.h) || !std::isfinite(s.gamma) ||
      !std::isfinite(s.xcg) || s.xcg != model.xcg) {
    return tr;
  }
  auto state = [&](const Vec3& z) {
    Vec<13> x = Vec<13>::Zero();
    x(0) = s.V;
    x(1) = z(0);
    x(4) = z(0) + s.gamma;
    x(11) = s.h;
    const F16Input u{z(1), z(2), 0, 0};
    x(12) = ideal.evaluate(0, x17_from_x13(ideal, x, u), u).power_cmd;
    return x;
  };
  auto residual = [&](const Vec3& z) -> Vec3 {
    const auto& d = ideal.data;
    if (z(0) * d.rtod * units::deg < d.cz.x().front() || z(0) * d.rtod * units::deg > d.cz.x().back() ||
        z(1) < d.throttle_min || z(1) > d.throttle_max || z(2) < d.elevator.pos_min ||
        z(2) > d.elevator.pos_max) {
      throw std::domain_error("trim trial outside alpha/throttle/elevator envelope");
    }
    const Vec<13> dx = f13(ideal, state(z), {z(1), z(2), 0, 0});
    return Vec3(dx(0) / ideal.grav.gravity(s.h), dx(1), dx(7));
  };
  try {
    const auto sol = solve_trim_system(residual, Vec3(5 * units::deg, 0.2, -units::deg));
    tr.ok = sol.ok;
    tr.iterations = sol.iterations;
    tr.residual = sol.residual;
    tr.alpha = sol.z(0);
    tr.throttle = sol.z(1);
    tr.de = sol.z(2);
    tr.u = {tr.throttle, tr.de, 0, 0};
    tr.x = x17_from_x13(model, state(sol.z), tr.u);
    if (tr.ok) {
      F16Model checked = model;
      LookupMonitor monitor;
      checked.monitor = &monitor;
      checked.evaluate(0, tr.x, tr.u);
      tr.ok = !monitor.out_of_domain;
      tr.message = tr.ok ? "converged" : monitor.first_violation;
    }
  } catch (const std::exception& e) {
    tr.message = e.what();
  }
  return tr;
}
}  // namespace flightsim
