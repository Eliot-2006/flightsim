#pragma once

#include <Eigen/Dense>
#include <array>
#include <string>
#include <vector>

namespace flightsim {

// Optional diagnostics for table lookups. Not thread-safe; pass by pointer, null is allowed.
struct LookupMonitor {
  bool out_of_domain = false;
  std::string
      first_violation;  // "<table>: <axis>=<value> outside [lo,hi]" for the first violation since reset()
  bool record_cells = false;  // set true to log every lookup
  // One entry per lookup when record_cells: {table-id hash, cell index}. Cell index is i for 1D and
  // i*(ny-1)+j for 2D, where i/j are the lower breakpoint indices of the (boundary-clamped) cell used.
  std::vector<std::array<int, 2>> cells;
  void reset();  // clears out_of_domain, first_violation and cells; keeps record_cells
};

// Semantics (Stevens & Lewis): piecewise linear / bilinear. Outside the breakpoint range the boundary cell is
// extended linearly (no distance limit) and the monitor, if given, flags out_of_domain. NaN inputs are
// flagged too. Breakpoints must be strictly increasing with at least two points (checked in the constructor,
// throws std::invalid_argument).
class Table1D {
 public:
  Table1D(std::string name, std::vector<double> x, std::vector<double> y);
  double operator()(double x, LookupMonitor* mon = nullptr) const;
  const std::string& name() const { return name_; }
  int id() const { return id_; }  // matches LookupMonitor::cells, for diagnostic classification
  const std::vector<double>& x() const { return x_; }

 private:
  std::string name_;
  int id_;
  std::vector<double> x_, y_;
};

// v(i,j) is the value at (x_i, y_j): rows follow x, columns follow y.
class Table2D {
 public:
  Table2D(std::string name, std::vector<double> x, std::vector<double> y, Eigen::MatrixXd v);
  double operator()(double x, double y, LookupMonitor* mon = nullptr) const;
  const std::string& name() const { return name_; }
  int id() const { return id_; }  // matches LookupMonitor::cells, for diagnostic classification
  const std::vector<double>& x() const { return x_; }
  const std::vector<double>& y() const { return y_; }

 private:
  std::string name_;
  int id_;
  std::vector<double> x_, y_;
  Eigen::MatrixXd v_;
};

}  // namespace flightsim
