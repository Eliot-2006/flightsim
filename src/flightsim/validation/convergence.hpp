#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "flightsim/core/integrators.hpp"
#include "flightsim/validation/cases.hpp"

namespace flightsim::validation {

struct ErrorNorms {
  double pos, vel, att, omega;  // Euclidean norms; att is the sign-invariant angle (PRD 4.2), rad
};
ErrorNorms state_error(const Vec<13>& x, const Vec<13>& ref);

struct RunResult {
  Vec<13> x_final;
  long nfev;
  double max_energy_drift = 0, max_momentum_drift = 0, max_qnorm_drift = 0;
};
// Fixed-step integration to c.t_final; h must be positive, finite and divide t_final.
RunResult integrate_case(const CaseDef& c, Integrator integ, double h, bool renormalize = true);

struct ConvergenceRow {
  std::string case_name;
  Integrator integ;
  double h;
  ErrorNorms err;
  long nfev;
  double err_E = 0, err_H = 0, err_qnorm = 0;  // Relative invariant drifts, absolute quaternion norm drift.
};
// h_k = h0 / 2^k. Uses an analytic reference, or independent RK4 refined per PRD 6.3.
// Throws if reference accuracy cannot be established (including roundoff-limited studies).
std::vector<ConvergenceRow> convergence_study(const CaseDef& c, Integrator integ, double h0, int levels,
                                              bool renormalize = true);

// Parse a study TOML and write convergence.csv, orders.json and norm_conservation.csv.
void write_convergence_study(const std::filesystem::path& study, const std::filesystem::path& out);

// Least-squares slope of log(err) vs log(h) over the asymptotic window (PRD 6.3): the finest run of levels
// whose error exceeds 1e3*eps*scale (above the roundoff floor) and whose successive local orders p_k agree
// within 10%. order is NaN (levels -1) when fewer than two usable levels exist.
struct OrderFit {
  double order;
  int first_level, last_level;
};
OrderFit observed_order(const std::vector<double>& h, const std::vector<double>& err, double scale);

}  // namespace flightsim::validation
