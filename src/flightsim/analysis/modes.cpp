#include "flightsim/analysis/modes.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <stdexcept>

namespace flightsim {
std::vector<Mode> identify_modes(const LinearModel& lm) {
  auto roots = [&](const std::array<int, 4>& ix) {
    Eigen::Matrix4d a;
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        a(i, j) = lm.A(ix[static_cast<std::size_t>(i)], ix[static_cast<std::size_t>(j)]);
      }
    }
    Eigen::EigenSolver<Eigen::Matrix4d> eig(a, false);
    if (eig.info() != Eigen::Success) {
      throw std::runtime_error("mode eigensolver failed");
    }
    return eig.eigenvalues().eval();
  };
  std::vector<std::complex<double>> lon, lat, real;
  // Power is held at trim: the secondary reference's longitudinal model has no power state.
  for (auto e : roots({0, 1, 4, 7})) {
    if (e.imag() >= 0) {
      lon.push_back(e);
    }
  }
  for (auto e : roots({2, 3, 6, 8})) {
    if (e.imag() > 0) {
      lat.push_back(e);
    }
    // Eigen returns exactly zero imaginary part for real roots of this real matrix.
    if (e.imag() == 0) {
      real.push_back(e);
    }
  }
  if (lon.empty() || lat.size() != 1 || real.size() != 2) {
    throw std::runtime_error("mode topology has no longitudinal roots or unsupported lateral roots");
  }
  std::sort(lon.begin(), lon.end(), [](auto a, auto b) { return std::abs(a) < std::abs(b); });
  std::sort(real.begin(), real.end(), [](auto a, auto b) { return std::abs(a) < std::abs(b); });
  std::vector<Mode> modes;
  auto add = [&](const char* name, std::complex<double> e) {
    modes.push_back({name, e, std::abs(e), -e.real() / std::abs(e), -1 / e.real()});
  };
  int aperiodic = 0;
  for (const auto& e : lon) {
    if (e.imag() > 0) {
      add(e == lon.front() ? "phugoid" : "short_period", e);
    } else {
      // At aft CG the short-period pair splits into real roots; retain the unstable root explicitly.
      const std::string name = "longitudinal_real_" + std::to_string(++aperiodic);
      add(name.c_str(), e);
    }
  }
  add("dutch_roll", lat[0]);
  add("roll", real[1]);
  add("spiral", real[0]);
  return modes;
}
}  // namespace flightsim
