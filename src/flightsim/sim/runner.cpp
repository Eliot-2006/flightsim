#include "flightsim/sim/runner.hpp"

#include <cmath>
#include <cstdio>

namespace flightsim::sim {

const char* to_string(StopReason r) {
  switch (r) {
    case StopReason::FinalTime:
      return "final_time";
    case StopReason::GroundContact:
      return "ground_contact";
    case StopReason::NonFinite:
      return "non_finite";
  }
  return "?";
}

namespace {

std::string describe_nonfinite(double t, long step, const Vec<13>& x) {
  std::string s = "non-finite state at t=";
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", t);
  s += buf;
  s += " s (step " + std::to_string(step) + "): x = [";
  for (int i = 0; i < 13; ++i) {
    std::snprintf(buf, sizeof buf, "%.17g", x(i));
    s += (i ? ", " : "");
    s += buf;
  }
  return s + "]";
}

}  // namespace

RunOutput run(const RigidBodyScenario& s) {
  const long n = std::lround(s.run.t_final / s.run.dt);
  const auto time_at = [&](long k) { return k == n ? s.run.t_final : static_cast<double>(k) * s.run.dt; };
  auto f = [&s](double t, const Vec<13>& x) {
    ForceMoment fm;
    for (const auto& force : s.forces) fm += force(t, x);
    return rigid_body_deriv(x, s.mp, fm, s.g_N);
  };

  RunOutput out;
  out.samples.reserve(static_cast<std::size_t>(n / s.run.output_every) + 2);
  Vec<13> x = s.x0;
  StepStats stats;
  out.samples.push_back({0.0, x});
  out.t_stop = s.run.t_final;
  for (long k = 0; k < n; ++k) {
    x = step(s.run.integ, f, time_at(k), x, s.run.dt, s.run.renormalize, stats);
    const double t = time_at(k + 1);
    if (!x.allFinite()) {
      out.reason = StopReason::NonFinite;
      out.t_stop = t;
      out.message = describe_nonfinite(t, k + 1, x);
      break;
    }
    if (s.run.stop_at_ground && -x(rb::kPos + 2) <= 0.0) {
      out.reason = StopReason::GroundContact;
      out.t_stop = t;
      out.samples.push_back({t, x});
      break;
    }
    if ((k + 1) % s.run.output_every == 0) out.samples.push_back({t, x});
  }
  out.nfev = stats.nfev;
  return out;
}

}  // namespace flightsim::sim
