#pragma once
#include <complex>
#include <string>
#include <vector>

#include "flightsim/analysis/linearize.hpp"
namespace flightsim {
struct Mode {
  std::string name;
  std::complex<double> eig;
  double wn, zeta, tau;
};
std::vector<Mode> identify_modes(const LinearModel&);
}  // namespace flightsim
