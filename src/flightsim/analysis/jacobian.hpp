#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "flightsim/math/quat.hpp"

namespace flightsim {

// Scaled central differences; scale_i is a positive characteristic size for state component i.
template <int N, int M, class F>
Eigen::Matrix<double, N, M> jacobian_central(F&& f, const Vec<M>& x, const Vec<M>& scale, double eps = 6e-6) {
  if (!x.allFinite() || !scale.allFinite() || (scale.array() <= 0).any() || !std::isfinite(eps) || eps <= 0)
    throw std::invalid_argument("jacobian_central requires finite x and positive finite scales/eps");
  Eigen::Matrix<double, N, M> J;
  for (int i = 0; i < M; ++i) {
    const double delta = eps * std::max(std::abs(x(i)), scale(i));
    Vec<M> plus = x, minus = x;
    plus(i) += delta;
    minus(i) -= delta;
    if (!plus.allFinite() || !minus.allFinite() || plus(i) == minus(i))
      throw std::invalid_argument("jacobian_central stencil is not representable");
    const Vec<N> fp = f(plus), fm = f(minus);
    if (!fp.allFinite() || !fm.allFinite()) throw std::runtime_error("jacobian_central non-finite f");
    J.col(i) = (fp - fm) / (2 * delta);
  }
  return J;
}
}  // namespace flightsim
