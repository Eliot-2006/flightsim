#pragma once

#include "flightsim/aircraft/f16_model.hpp"

namespace flightsim {

// S&L chart, SI: [V, alpha, beta, phi, theta, psi, p, q, r, pN, pE, h, P].
// V is AIR-relative, altitude is up. Surface states initialize from u, in rad.
Vec<17> x17_from_x13(const F16Model&, const Vec<13>& x13, const F16Input& u, double t = 0);
Vec<13> x13_from_x17(const F16Model&, const Vec<17>& x, double t = 0);
// Assumes steady wind (no wind time derivative). Chart singular at Euler pitch/sideslip poles.
// At V < 1e-9, V/alpha/beta derivatives are defined as zero; use x17 dynamics to depart from rest.
Vec<13> f13(const F16Model&, const Vec<13>& x13, const F16Input& u, double t = 0);

}  // namespace flightsim
