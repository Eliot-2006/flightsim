#pragma once

namespace flightsim {
struct JacobiSnCnDn {
  double sn, cn, dn;
};
// Jacobi elliptic functions, parameter m = k^2, finite u and 0 <= m < 1 (A&S 16.4).
JacobiSnCnDn jacobi_elliptic(double u, double m);
}  // namespace flightsim
