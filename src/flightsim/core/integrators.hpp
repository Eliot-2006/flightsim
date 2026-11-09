#pragma once

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "flightsim/core/state.hpp"

namespace flightsim {

enum class Integrator { Euler, RK4, SymplecticEuler, LieEuler, RKMK4 };

inline const char* to_string(Integrator i) {
  switch (i) {
    case Integrator::Euler:
      return "euler";
    case Integrator::RK4:
      return "rk4";
    case Integrator::SymplecticEuler:
      return "symplectic_euler";
    case Integrator::LieEuler:
      return "lie_euler";
    case Integrator::RKMK4:
      return "rkmk4";
  }
  return "?";
}

inline Integrator integrator_from_string(std::string_view s) {
  for (Integrator i : {Integrator::Euler, Integrator::RK4, Integrator::SymplecticEuler, Integrator::LieEuler,
                       Integrator::RKMK4})
    if (s == to_string(i)) return i;
  throw std::invalid_argument("unknown integrator '" + std::string(s) +
                              "' (expected euler | rk4 | symplectic_euler | lie_euler | rkmk4)");
}

struct StepStats {
  long nfev = 0;  // derivative evaluations (work-precision x-axis)
};

// f: (double t, const Vec<N>& x) -> Vec<N>. Fixed-size, no heap allocation.
template <int N, class F>
Vec<N> euler_step(F&& f, double t, const Vec<N>& x, double h) {
  const Vec<N> k1 = f(t, x);
  return x + h * k1;
}

template <int N, class F>
Vec<N> rk4_step(F&& f, double t, const Vec<N>& x, double h) {
  const Vec<N> k1 = f(t, x);
  const Vec<N> k2 = f(t + 0.5 * h, Vec<N>(x + 0.5 * h * k1));
  const Vec<N> k3 = f(t + 0.5 * h, Vec<N>(x + 0.5 * h * k2));
  const Vec<N> k4 = f(t + h, Vec<N>(x + h * k3));
  return x + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

namespace detail {
inline Vec4 lie_attitude_update(const Vec4& q, const Vec3& theta) {
  const Vec4 e = quat_from_rotvec(theta);
  if (e(0) < 0) return quat_mul(q, e);
  // q + q*(exp(theta/2)-1) avoids a rounded cosine bias on every small step; no state projection.
  Vec4 delta = e;
  delta(0) = -e.tail<3>().squaredNorm() / (1.0 + e(0));
  return q + quat_mul(q, delta);
}
}  // namespace detail

// Semi-implicit only for separable translation; rotation and extra states use explicit Euler.
template <int N, class F>
Vec<N> symplectic_euler_step(F&& f, double t, const Vec<N>& x, double h) {
  static_assert(N >= rb::kN);
  Vec<N> y = euler_step(f, t, x, h);
  rb::pos(y) = rb::pos(x) + h * rb::vel(y);
  return y;
}

template <int N, class F>
Vec<N> lie_euler_step(F&& f, double t, const Vec<N>& x, double h) {
  static_assert(N >= rb::kN);
  Vec<N> y = euler_step(f, t, x, h);
  rb::quat(y) = detail::lie_attitude_update(rb::quat(x), h * rb::omega(x));
  return y;
}

template <int N, class F>
Vec<N> rkmk4_step(F&& f, double t, const Vec<N>& x, double h) {
  static_assert(N >= rb::kN);
  auto stage = [&](double time, Vec<N> y, const Vec3& theta) {
    rb::quat(y) = detail::lie_attitude_update(rb::quat(x), theta);
    const Vec3 w = rb::omega(y);
    Vec<N> d = f(time, y);
    rb::quat(d).setZero();
    // Body-frame right multiplication requires dexp^{-1}_{-theta}, hence +1/2 (PRD 6.1).
    return std::pair<Vec<N>, Vec3>{d, w + 0.5 * theta.cross(w) + theta.cross(theta.cross(w)) / 12.0};
  };
  const auto [k1, a1] = stage(t, x, Vec3::Zero());
  const auto [k2, a2] = stage(t + h / 2, x + h / 2 * k1, h / 2 * a1);
  const auto [k3, a3] = stage(t + h / 2, x + h / 2 * k2, h / 2 * a2);
  const auto [k4, a4] = stage(t + h, x + h * k3, h * a3);
  Vec<N> y = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  rb::quat(y) = detail::lie_attitude_update(rb::quat(x), h / 6 * (a1 + 2 * a2 + 2 * a3 + a4));
  return y;
}

// Projects the quaternion block back to unit norm (PRD 4.2 policy for integrators that step q in R^4).
template <int N>
void renormalize_quat(Vec<N>& x) {
  static_assert(N >= rb::kN);
  x.template segment<4>(rb::kQuat).normalize();
}

// One step on a rigid-body-layout state. Lie methods never project the state quaternion.
template <int N, class F>
Vec<N> step(Integrator integ, F&& f, double t, const Vec<N>& x, double h, bool renormalize,
            StepStats& stats) {
  static_assert(N >= rb::kN);
  Vec<N> y;
  switch (integ) {
    case Integrator::Euler:
      y = euler_step(f, t, x, h);
      stats.nfev += 1;
      break;
    case Integrator::RK4:
      y = rk4_step(f, t, x, h);
      stats.nfev += 4;
      break;
    case Integrator::SymplecticEuler:
      y = symplectic_euler_step(f, t, x, h);
      stats.nfev += 1;
      break;
    case Integrator::LieEuler:
      y = lie_euler_step(f, t, x, h);
      stats.nfev += 1;
      return y;
    case Integrator::RKMK4:
      y = rkmk4_step(f, t, x, h);
      stats.nfev += 4;
      return y;
    default:
      throw std::invalid_argument("invalid integrator");
  }
  if (renormalize) renormalize_quat(y);
  return y;
}

}  // namespace flightsim
