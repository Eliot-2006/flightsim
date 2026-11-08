#pragma once

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
};
// Fixed-step integration to c.t_final; h must divide t_final. Throws for integrators not implemented yet.
RunResult integrate_case(const CaseDef& c, Integrator integ, double h, bool renormalize = true);

struct ConvergenceRow {
  std::string case_name;
  Integrator integ;
  double h;
  ErrorNorms err;
  long nfev;
};
// h_k = h0 / 2^k, k = 0..levels-1, error at t_final against c.exact (throws if the case has no closed form).
std::vector<ConvergenceRow> convergence_study(const CaseDef& c, Integrator integ, double h0, int levels);

// Least-squares slope of log(err) vs log(h) over the asymptotic window (PRD 6.3): the finest run of levels
// whose error exceeds 1e3*eps*scale (above the roundoff floor) and whose successive local orders p_k agree
// within 10%. order is NaN (levels -1) when fewer than two usable levels exist.
struct OrderFit {
  double order;
  int first_level, last_level;
};
OrderFit observed_order(const std::vector<double>& h, const std::vector<double>& err, double scale);

}  // namespace flightsim::validation
