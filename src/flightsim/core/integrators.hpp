#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

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

// Projects the quaternion block back to unit norm (PRD 4.2 policy for integrators that step q in R^4).
template <int N>
void renormalize_quat(Vec<N>& x) {
  static_assert(N >= rb::kN);
  x.template segment<4>(rb::kQuat).normalize();
}

// One step of the chosen method on a rigid-body-layout state. Integrators beyond Euler/RK4 arrive in M2.
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
    default:
      throw std::logic_error(std::string("integrator '") + to_string(integ) + "' is not implemented yet");
  }
  if (renormalize) renormalize_quat(y);
  return y;
}

}  // namespace flightsim
