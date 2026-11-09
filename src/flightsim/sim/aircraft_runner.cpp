#include "flightsim/sim/aircraft_runner.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/version.hpp"

#ifndef FLIGHTSIM_VERSION
#error "FLIGHTSIM_VERSION must be defined by the build"
#endif

namespace flightsim::sim {
namespace {
F16Input input_at(const AircraftScenario& s, const F16Input& base, double t) {
  F16Input u = base;
  for (const auto& e : s.inputs) {
    if (e.t_start_s <= t && t < e.t_end_s) {
      (e.channel == InputEvent::Throttle ? u.throttle
       : e.channel == InputEvent::De     ? u.de
       : e.channel == InputEvent::Da     ? u.da
                                         : u.dr) += e.value;
    }
  }
  return u;
}
}  // namespace

AircraftRunOutput run_aircraft(const AircraftScenario& s) {
  if (!s.run) {
    throw std::runtime_error("aircraft scenario has no [sim] table");
  }
  const auto& rc = *s.run;
  AircraftRunOutput out;
  Vec<17> x;
  if (s.from_trim) {
    const auto tr = trim_level(s.model, s.trim);
    if (!tr.ok) {
      throw std::runtime_error("trim failed: " + tr.message);
    }
    x = tr.x;
    out.u_trim = tr.u;
  } else {
    out.u_trim = s.u_base;
    x = x17_from_x13(s.model, s.x13_0, s.u_base);
  }
  LookupMonitor monitor;
  F16Model model = s.model;
  model.monitor = &monitor;

  const long n = std::lround(rc.t_final / rc.dt);
  const auto time_at = [&](long k) { return k == n ? rc.t_final : static_cast<double>(k) * rc.dt; };
  StepStats stats;
  const auto record = [&](double t, const Vec<17>& xs, const F16Input& u) {
    out.t.push_back(t);
    out.x.push_back(xs);
    out.u.push_back(u);
  };
  F16Input u = input_at(s, out.u_trim, 0.0);
  record(0.0, x, u);
  out.t_stop = rc.t_final;
  for (long k = 0; k < n; ++k) {
    const double t = time_at(k);
    u = input_at(s, out.u_trim, t);
    // M7 inserts the fixed-rate controller call here: it may overwrite `u` before the step.
    try {
      x = step(
          rc.integ, [&](double tt, const Vec<17>& xx) { return model.deriv(tt, xx, u); }, t, x, rc.dt, true, stats);
    } catch (const std::exception& e) {  // atmosphere altitude domain, non-finite inputs
      out.reason = StopReason::NonFinite;
      out.t_stop = time_at(k + 1);
      out.message = std::string("step failed at t=") + std::to_string(t) + " s: " + e.what();
      break;
    }
    const double tn = time_at(k + 1);
    if (!x.allFinite()) {
      out.reason = StopReason::NonFinite;
      out.t_stop = tn;
      out.message = "non-finite state at t=" + std::to_string(tn) + " s";
      break;
    }
    const bool ground = rc.stop_at_ground && -x(rb::kPos + 2) <= 0.0;
    if (ground || (k + 1) % rc.output_every == 0) {
      record(tn, x, input_at(s, out.u_trim, tn));
    }
    if (ground) {
      out.reason = StopReason::GroundContact;
      out.t_stop = tn;
      break;
    }
  }
  out.out_of_domain = monitor.out_of_domain;
  out.first_violation = monitor.first_violation;
  return out;
}

namespace {
constexpr const char* kCols[] = {
    "t_s", "pN_m", "pE_m", "pD_m", "vN_mps", "vE_mps", "vD_mps", "q0", "q1", "q2", "q3", "p_radps", "q_radps", "r_radps",
    "power_pct", "de_state_rad", "da_state_rad", "dr_state_rad", "V_mps", "alpha_rad", "beta_rad", "phi_rad",
    "theta_rad", "psi_rad", "h_m", "mach", "qbar_Pa", "thrust_N", "throttle", "de_rad", "da_rad", "dr_rad",
    "V_fps", "alpha_deg", "beta_deg", "phi_deg", "theta_deg", "psi_deg", "p_dps", "q_dps", "r_dps", "h_ft",
    "pN_ft", "pE_ft", "de_deg", "da_deg", "dr_deg"};
}  // namespace

void write_aircraft_csv(const std::filesystem::path& path, const AircraftScenario& s, const AircraftRunOutput& out) {
  std::ofstream f(path);
  if (!f) {
    throw std::runtime_error("cannot write " + path.string());
  }
  F16Model model = s.model;
  model.monitor = nullptr;
  bool first = true;
  for (const char* c : kCols) {
    f << (first ? "" : ",") << c;
    first = false;
  }
  f << '\n';
  char buf[40];
  for (std::size_t i = 0; i < out.t.size(); ++i) {
    const Vec<17>& x = out.x[i];
    const F16Input& u = out.u[i];
    const auto ev = model.evaluate(out.t[i], x, u);
    const Vec<13> c = x13_from_x17(model, x, out.t[i]);
    const double d = units::deg, ft = units::ft;
    const double row[] = {out.t[i], x(0), x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), x(9), x(10), x(11), x(12),
                          x(13), x(14), x(15), x(16), c(0), c(1), c(2), c(3), c(4), c(5), c(11), ev.mach, ev.qbar,
                          ev.thrust, u.throttle, u.de, u.da, u.dr, c(0) / ft, c(1) / d, c(2) / d, c(3) / d,
                          c(4) / d, c(5) / d, c(6) / d, c(7) / d, c(8) / d, c(11) / ft, c(9) / ft, c(10) / ft,
                          u.de / d, u.da / d, u.dr / d};
    for (std::size_t j = 0; j < std::size(row); ++j) {
      std::snprintf(buf, sizeof buf, "%.17g", row[j]);
      f << (j ? "," : "") << buf;
    }
    f << '\n';
  }
}

void write_aircraft_metadata(const std::filesystem::path& path, const AircraftScenario& s,
                             const AircraftRunOutput& out) {
  const auto& rc = *s.run;
  toml::table m;
  m.insert("version", FLIGHTSIM_VERSION);
  m.insert("scenario", s.source_path);
  m.insert("integrator", to_string(rc.integ));
  m.insert("dt_s", rc.dt);
  m.insert("t_final_s", rc.t_final);
  m.insert("xcg", s.model.xcg);
  m.insert("stop_reason", to_string(out.reason));
  m.insert("t_stop_s", out.t_stop);
  m.insert("out_of_domain", out.out_of_domain);
  m.insert("first_violation", out.first_violation);
  m.insert("message", out.message);
  std::ofstream f(path);
  if (!f) {
    throw std::runtime_error("cannot write " + path.string());
  }
  f << toml::json_formatter{m} << '\n';
}
}  // namespace flightsim::sim
