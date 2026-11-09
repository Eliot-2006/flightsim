#include "flightsim/math/interp.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>

namespace flightsim {

void LookupMonitor::reset() {
  out_of_domain = false;
  first_violation.clear();
  cells.clear();
}

namespace {

int table_id(const std::string& name) {
  return static_cast<int>(std::hash<std::string>{}(name) & 0x7fffffffu);
}

void check_axis(const std::string& name, const char* axis, const std::vector<double>& x) {
  if (x.size() < 2) throw std::invalid_argument(name + ": axis " + axis + " needs at least 2 breakpoints");
  for (std::size_t i = 1; i < x.size(); ++i)
    if (!(x[i] > x[i - 1]))  // also rejects NaN
      throw std::invalid_argument(name + ": axis " + axis + " breakpoints not strictly increasing at index " +
                                  std::to_string(i));
}

// Lower breakpoint index of the cell used for v: boundary cell when v is outside (or NaN).
std::size_t find_cell(const std::vector<double>& x, double v) {
  const auto it = std::upper_bound(x.begin(), x.end(), v);
  const auto idx = static_cast<std::size_t>(it - x.begin());
  return std::min(idx == 0 ? std::size_t{0} : idx - 1, x.size() - 2);
}

// Returns the fractional position in cell i (outside [0,1] when extrapolating) and flags the monitor.
double cell_fraction(const std::string& name, const char* axis, const std::vector<double>& x, std::size_t i,
                     double v, LookupMonitor* mon) {
  if (mon && !mon->out_of_domain && !(v >= x.front() && v <= x.back())) {
    std::ostringstream os;
    os.precision(10);
    os << name << ": " << axis << "=" << v << " outside [" << x.front() << "," << x.back() << "]";
    mon->out_of_domain = true;
    mon->first_violation = os.str();
  }
  return (v - x[i]) / (x[i + 1] - x[i]);
}

}  // namespace

Table1D::Table1D(std::string name, std::vector<double> x, std::vector<double> y)
    : name_(std::move(name)), id_(table_id(name_)), x_(std::move(x)), y_(std::move(y)) {
  check_axis(name_, "x", x_);
  if (y_.size() != x_.size()) throw std::invalid_argument(name_ + ": y size does not match x size");
}

double Table1D::operator()(double x, LookupMonitor* mon) const {
  const std::size_t i = find_cell(x_, x);
  const double t = cell_fraction(name_, "x", x_, i, x, mon);
  if (mon && mon->record_cells) mon->cells.push_back({id_, static_cast<int>(i)});
  // (1-t)a + t b is exact at t = 0 and t = 1, so breakpoints return the stored value bit-for-bit.
  return (1.0 - t) * y_[i] + t * y_[i + 1];
}

Table2D::Table2D(std::string name, std::vector<double> x, std::vector<double> y, Eigen::MatrixXd v)
    : name_(std::move(name)), id_(table_id(name_)), x_(std::move(x)), y_(std::move(y)), v_(std::move(v)) {
  check_axis(name_, "x", x_);
  check_axis(name_, "y", y_);
  if (static_cast<std::size_t>(v_.rows()) != x_.size() || static_cast<std::size_t>(v_.cols()) != y_.size())
    throw std::invalid_argument(name_ +
                                ": value matrix shape does not match breakpoints (x rows, y columns)");
}

double Table2D::operator()(double x, double y, LookupMonitor* mon) const {
  const std::size_t i = find_cell(x_, x);
  const std::size_t j = find_cell(y_, y);
  const double s = cell_fraction(name_, "x", x_, i, x, mon);
  const double t = cell_fraction(name_, "y", y_, j, y, mon);
  if (mon && mon->record_cells) mon->cells.push_back({id_, static_cast<int>(i * (y_.size() - 1) + j)});
  const auto ii = static_cast<Eigen::Index>(i);
  const auto jj = static_cast<Eigen::Index>(j);
  const double a = (1.0 - t) * v_(ii, jj) + t * v_(ii, jj + 1);
  const double b = (1.0 - t) * v_(ii + 1, jj) + t * v_(ii + 1, jj + 1);
  return (1.0 - s) * a + s * b;
}

}  // namespace flightsim
