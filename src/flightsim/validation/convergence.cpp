#include "flightsim/validation/convergence.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <toml++/toml.hpp>

#include "flightsim/math/quat.hpp"
#include "flightsim/validation/analytic.hpp"

namespace flightsim::validation {

ErrorNorms state_error(const Vec<13>& x, const Vec<13>& ref) {
  return {(rb::pos(x) - rb::pos(ref)).norm(), (rb::vel(x) - rb::vel(ref)).norm(),
          attitude_error(Vec4(rb::quat(ref)), Vec4(rb::quat(x))), (rb::omega(x) - rb::omega(ref)).norm()};
}

RunResult integrate_case(const CaseDef& c, Integrator integ, double h, bool renormalize) {
  if (!std::isfinite(h) || h <= 0 || !std::isfinite(c.t_final) || c.t_final <= 0 ||
      c.t_final / h >= static_cast<double>(std::numeric_limits<long>::max()))
    throw std::invalid_argument(c.name + ": invalid integration time grid");
  const long n = std::lround(c.t_final / h);
  if (n < 1 || std::abs(static_cast<double>(n) * h - c.t_final) > 1e-9 * c.t_final)
    throw std::invalid_argument(c.name + ": step h does not divide t_final");
  auto f = [&c](double t, const Vec<13>& x) {
    ForceMoment fm;
    for (const auto& force : c.forces) fm += force(t, x);
    return rigid_body_deriv(x, c.mp, fm, c.g_N);
  };
  Vec<13> x = c.x0;
  StepStats stats;
  for (long i = 0; i < n; ++i) {
    x = step(integ, f, static_cast<double>(i) * h, x, h, renormalize, stats);
    if (!x.allFinite()) throw std::runtime_error(c.name + ": non-finite state during integration");
  }
  return {x, stats.nfev};
}

namespace {
bool reference_accurate(const ErrorNorms& difference, const ErrorNorms& error) {
  // RK4 step-halving estimates fine-grid error as |x_h-x_2h|/(2^4-1).
  auto accurate = [](double d, double e) { return d == 0 || d / 15 < 1e-3 * e; };
  return accurate(difference.pos, error.pos) && accurate(difference.vel, error.vel) &&
         accurate(difference.att, error.att) && accurate(difference.omega, error.omega);
}
}  // namespace

std::vector<ConvergenceRow> convergence_study(const CaseDef& c, Integrator integ, double h0, int levels,
                                              bool renormalize) {
  if (levels < 2 || levels > 30) throw std::invalid_argument("convergence levels must be in [2,30]");
  std::vector<RunResult> runs;
  for (int k = 0; k < levels; ++k) runs.push_back(integrate_case(c, integ, std::ldexp(h0, -k), renormalize));
  Vec<13> ref;
  if (c.exact)
    ref = c.exact(c.t_final);
  else {
    double href = std::ldexp(h0, -(levels - 1)) / 16;
    Vec<13> coarse = integrate_case(c, Integrator::RK4, 2 * href).x_final;
    ref = integrate_case(c, Integrator::RK4, href).x_final;
    int refinements = 0;
    for (;;) {
      ErrorNorms estimate = state_error(ref, coarse);
      ErrorNorms finest = state_error(runs.back().x_final, ref);
      if (c.exact_omega) {
        estimate.omega = 0;
        finest.omega = 0;
      }
      if (reference_accurate(estimate, finest)) break;
      if (++refinements > 8)
        throw std::runtime_error(
            c.name + ": RK4 fine-reference accuracy not established; use a coarser study above roundoff");
      coarse = ref;
      href /= 2;
      ref = integrate_case(c, Integrator::RK4, href).x_final;
    }
  }
  if (c.exact_omega) rb::omega(ref) = c.exact_omega(c.t_final);
  auto energy = [&](const Vec<13>& x) { return c.energy ? c.energy(x) : rotational_energy(x, c.mp); };
  const double E0 = energy(c.x0);
  const Vec3 H0 = angular_momentum_N(c.x0, c.mp);
  std::vector<ConvergenceRow> rows;
  for (int k = 0; k < levels; ++k) {
    const auto& r = runs[static_cast<std::size_t>(k)];
    rows.push_back({c.name, integ, std::ldexp(h0, -k), state_error(r.x_final, ref), r.nfev,
                    E0 > 0 ? std::abs(energy(r.x_final) / E0 - 1) : 0,
                    H0.norm() > 0 ? (angular_momentum_N(r.x_final, c.mp) - H0).norm() / H0.norm() : 0,
                    std::abs(rb::quat(r.x_final).norm() - rb::quat(c.x0).norm())});
  }
  return rows;
}

namespace {
CaseDef named_case(const std::string& name) {
  for (auto factory :
       {case_T0, case_T1, case_T2, case_T3, case_T4, case_T5, case_R1, case_R2, case_R3, case_R4, case_R5}) {
    auto c = factory();
    if (c.name == name) return c;
  }
  throw std::invalid_argument("unknown validation case '" + name + "'");
}

std::vector<std::string> string_array(const toml::table& t, const char* key) {
  const auto* a = t[key].as_array();
  if (!a || a->empty()) throw std::invalid_argument(std::string(key) + ": expected nonempty string array");
  std::vector<std::string> result;
  for (const auto& v : *a) {
    const auto s = v.value<std::string>();
    if (!s) throw std::invalid_argument(std::string(key) + ": expected string array");
    result.push_back(*s);
  }
  return result;
}

struct Study {
  std::vector<CaseDef> cases;
  std::vector<Integrator> integrators;
  double h0;
  int levels;
  bool renormalize;
};

std::vector<Study> load_studies(const std::filesystem::path& path) {
  const auto root = toml::parse_file(path.string());
  for (const auto& [key, value] : root)
    if (key != "studies") throw std::invalid_argument("unknown study key: " + std::string(key.str()));
  const auto* entries = root["studies"].as_array();
  if (!entries || entries->empty()) throw std::invalid_argument("expected nonempty [[studies]]");
  std::vector<Study> studies;
  std::set<std::pair<std::string, Integrator>> pairs;
  for (const auto& entry : *entries) {
    const auto* t = entry.as_table();
    if (!t) throw std::invalid_argument("studies: expected table");
    for (const auto& [key, value] : *t)
      if (key != "cases" && key != "integrators" && key != "h0" && key != "levels" && key != "renormalize")
        throw std::invalid_argument("unknown study key: " + std::string(key.str()));
    const auto h0 = (*t)["h0"].value<double>();
    const auto levels = (*t)["levels"].value<int>();
    const auto renormalize = (*t)["renormalize"].value<bool>();
    if (!h0 || !std::isfinite(*h0) || *h0 <= 0)
      throw std::invalid_argument("h0: expected positive finite number");
    if (!levels || *levels < 8 || *levels > 30)
      throw std::invalid_argument("levels: expected integer in [8,30]");
    if (t->contains("renormalize") && !renormalize)
      throw std::invalid_argument("renormalize: expected boolean");
    Study s{{}, {}, *h0, *levels, renormalize.value_or(true)};
    for (const auto& name : string_array(*t, "cases")) s.cases.push_back(named_case(name));
    for (const auto& name : string_array(*t, "integrators"))
      s.integrators.push_back(integrator_from_string(name));
    for (const auto& c : s.cases)
      for (auto method : s.integrators)
        if (!pairs.emplace(c.name, method).second)
          throw std::invalid_argument("duplicate case/integrator study");
    studies.push_back(std::move(s));
  }
  return studies;
}

std::vector<std::string> asserted_metrics(const std::string& name, Integrator method) {
  const bool ordinary = method == Integrator::Euler || method == Integrator::RK4;
  const bool lie = method == Integrator::LieEuler || method == Integrator::RKMK4;
  if ((ordinary && (name == "T2" || name == "T3" || name == "T4")) ||
      (method == Integrator::SymplecticEuler && name == "T4"))
    return {"pos", "vel"};
  if (name == "R1" && ordinary) return {method == Integrator::Euler ? "qnorm" : "att"};
  if (name == "R2" && (ordinary || lie)) return {"att", "omega"};
  if (name == "R3" && (ordinary || lie)) return {"omega", "H"};
  return {};
}

void write_text(const std::filesystem::path& path, const std::string& text) {
  std::ofstream stream;
  stream.exceptions(std::ios::failbit | std::ios::badbit);
  stream.open(path, std::ios::binary);
  stream << text;
  stream.close();
}
}  // namespace

void write_convergence_study(const std::filesystem::path& study_path, const std::filesystem::path& out) {
  const auto studies = load_studies(study_path);
  std::ostringstream csv, comparison;
  csv << std::setprecision(17)
      << "case,integrator,h,err_pos,err_vel,err_att,err_omega,err_E,err_H,nfev,err_qnorm,renormalize\n";
  comparison << std::setprecision(17) << "case,integrator,h,renormalize,err_qnorm,err_E,err_H,nfev\n";
  toml::array orders;
  bool all_passed = true;
  std::set<std::string> compared;
  for (const auto& study : studies)
    for (const auto& c : study.cases) {
      for (auto method : study.integrators) {
        const auto rows = convergence_study(c, method, study.h0, study.levels, study.renormalize);
        const auto asserted = asserted_metrics(c.name, method);
        if (c.name == "R1" && method == Integrator::Euler && study.renormalize)
          throw std::invalid_argument(
              "R1 Euler order 1 uses unprojected quaternion norm; set renormalize=false (angle "
              "superconverges at order 2)");
        std::array<std::vector<double>, 7> errors;
        const std::array<std::string, 7> names = {"pos", "vel", "att", "omega", "E", "H", "qnorm"};
        const Vec<13> scale_state = c.exact ? c.exact(c.t_final) : c.x0;
        const std::array<double, 7> scales = {std::max(1., rb::pos(scale_state).norm()),
                                              std::max(1., rb::vel(scale_state).norm()),
                                              1,
                                              std::max(1., rb::omega(scale_state).norm()),
                                              1,
                                              1,
                                              1};
        std::vector<double> h;
        for (const auto& row : rows) {
          h.push_back(row.h);
          const std::array<double, 7> e = {row.err.pos, row.err.vel, row.err.att,  row.err.omega,
                                           row.err_E,   row.err_H,   row.err_qnorm};
          for (std::size_t i = 0; i < e.size(); ++i) errors[i].push_back(e[i]);
          csv << c.name << ',' << to_string(method) << ',' << row.h << ',' << row.err.pos << ','
              << row.err.vel << ',' << row.err.att << ',' << row.err.omega << ',' << row.err_E << ','
              << row.err_H << ',' << row.nfev << ',' << row.err_qnorm << ',' << study.renormalize << '\n';
        }
        toml::array metrics;
        for (std::size_t i = 0; i < names.size(); ++i) {
          const auto fit = observed_order(h, errors[i], scales[i]);
          const bool required = std::find(asserted.begin(), asserted.end(), names[i]) != asserted.end();
          toml::table metric{{"metric", names[i]},
                             {"asserted", required},
                             {"first_level", fit.first_level},
                             {"last_level", fit.last_level}};
          if (std::isfinite(fit.order)) metric.insert("order", fit.order);
          if (required) {
            const bool fourth = method == Integrator::RK4 || method == Integrator::RKMK4;
            const double lower = fourth ? 3.9 : .95, upper = fourth ? 4.1 : 1.05;  // PRD 10.3 bounds.
            const bool pass = std::isfinite(fit.order) && fit.order >= lower && fit.order <= upper;
            metric.insert("lower_bound", lower);
            metric.insert("upper_bound", upper);
            metric.insert("within_bounds", pass);
            all_passed = all_passed && pass;
          }
          metrics.push_back(std::move(metric));
        }
        toml::table result{
            {"case", c.name}, {"integrator", to_string(method)}, {"renormalize", study.renormalize},
            {"h0", study.h0}, {"levels", study.levels},          {"metrics", std::move(metrics)}};
        if (c.name == "R1" && method == Integrator::Euler)
          result.insert(
              "note",
              "Order 1 is the unprojected quaternion norm error; constant-rate Euler angle is order 2.");
        orders.push_back(std::move(result));
      }
      if (c.name.starts_with("R") && compared.insert(c.name).second) {
        const double E0 = rotational_energy(c.x0, c.mp);
        const Vec3 H0 = angular_momentum_N(c.x0, c.mp);
        for (auto method : {Integrator::Euler, Integrator::RK4, Integrator::RKMK4}) {
          const auto r = integrate_case(c, method, .01, false);
          comparison << c.name << ',' << to_string(method) << ",0.01,false,"
                     << std::abs(rb::quat(r.x_final).norm() - 1) << ','
                     << std::abs(rotational_energy(r.x_final, c.mp) / E0 - 1) << ','
                     << (angular_momentum_N(r.x_final, c.mp) - H0).norm() / H0.norm() << ',' << r.nfev
                     << '\n';
        }
      }
    }
  toml::table result{{"all_assertions_passed", all_passed}, {"orders", std::move(orders)}};
  std::ostringstream json;
  json << toml::json_formatter{result} << '\n';
  std::filesystem::create_directories(out);
  write_text(out / "convergence.csv", csv.str());
  write_text(out / "orders.json", json.str());
  write_text(out / "norm_conservation.csv", comparison.str());
  if (!all_passed)
    throw std::runtime_error("convergence orders outside PRD 10.3 bounds; see " +
                             (out / "orders.json").string());
}

OrderFit observed_order(const std::vector<double>& h, const std::vector<double>& err, double scale) {
  const OrderFit none{std::numeric_limits<double>::quiet_NaN(), -1, -1};
  if (h.size() != err.size() || h.size() < 2) return none;
  const double floor_ = 1e3 * std::numeric_limits<double>::epsilon() * scale;

  // Finest level still above the roundoff floor (h decreases with the index).
  int last = -1;
  for (std::size_t i = 0; i < err.size(); ++i)
    if (err[i] > floor_) last = static_cast<int>(i);
  if (last < 1) return none;

  // Local orders p_k between level k and k+1; extend the window backwards from the finest level while they
  // agree.
  auto p = [&](int k) {
    const auto i = static_cast<std::size_t>(k);
    return std::log(err[i] / err[i + 1]) / std::log(h[i] / h[i + 1]);
  };
  int first = last - 1;
  while (first > 0 && err[static_cast<std::size_t>(first - 1)] > floor_ &&
         std::abs(p(first - 1) - p(first)) <= 0.1 * std::abs(p(first)))
    --first;

  // Least-squares slope of log(err) against log(h) over [first, last].
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  const double m = last - first + 1;
  for (int k = first; k <= last; ++k) {
    const double x = std::log(h[static_cast<std::size_t>(k)]), y = std::log(err[static_cast<std::size_t>(k)]);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
  }
  return {(m * sxy - sx * sy) / (m * sxx - sx * sx), first, last};
}

}  // namespace flightsim::validation
