#include "flightsim/math/elliptic.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace flightsim {

JacobiSnCnDn jacobi_elliptic(double u, double m) {
  if (!std::isfinite(u) || !std::isfinite(m) || m < 0 || m >= 1)
    throw std::invalid_argument("jacobi_elliptic requires finite u and 0 <= m < 1 (m=k^2)");
  if (m == 0) return {std::sin(u), std::cos(u), 1};

  // AGM followed by descending Landen transformations (A&S 16.4).
  std::array<double, 32> a{}, c{};
  a[0] = 1;
  double b = std::sqrt(1 - m);
  int n = 0;
  do {
    const auto i = static_cast<std::size_t>(n);
    a[i + 1] = (a[i] + b) / 2;
    c[i + 1] = (a[i] - b) / 2;
    b = std::sqrt(a[i] * b);
    ++n;
  } while (std::abs(c[static_cast<std::size_t>(n)]) >
           std::numeric_limits<double>::epsilon() * a[static_cast<std::size_t>(n)]);

  const double an = a[static_cast<std::size_t>(n)];
  const double period = 2 * std::numbers::pi / an;
  double phi = std::ldexp(an * std::remainder(u, period), n);
  for (int j = n; j > 0; --j) {
    const auto i = static_cast<std::size_t>(j);
    phi = (phi + std::asin(c[i] / a[i] * std::sin(phi))) / 2;
  }
  const double sn = std::sin(phi), cn = std::cos(phi);
  return {sn, cn, std::sqrt(1 - m * sn * sn)};
}
}  // namespace flightsim
