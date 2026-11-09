// PRD 11: F-16 benchmarks. Derivative evaluation, one RK4 step, the 60 s closed-loop C1 run, trim,
// linearization.
#include <benchmark/benchmark.h>

#include <filesystem>

#include "flightsim/analysis/linearize.hpp"
#include "flightsim/analysis/trim.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/sim/aircraft_runner.hpp"

namespace {
using namespace flightsim;
const std::filesystem::path root = FLIGHTSIM_SOURCE_DIR;

// The closed-loop operating point (scenarios/f16_ap_c1_alt100.toml): 502 ft/s, 10,000 ft, actuators and
// engine lag on.
struct Operating {
  sim::AircraftScenario scenario = sim::load_aircraft_scenario(root / "scenarios/f16_ap_c1_alt100.toml");
  TrimResult trim = trim_level(scenario.model, scenario.trim);
};
const Operating& operating() {
  static const Operating op;
  return op;
}

void BM_F16Deriv(benchmark::State& state) {
  const auto& op = operating();
  const double t = 0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(op.scenario.model.deriv(t, op.trim.x, op.trim.u));
  }
}
BENCHMARK(BM_F16Deriv);

// One RK4 step at the production dt, including the quaternion renormalization the runner applies. Every
// iteration starts from the trimmed state: stepping on would drift off trim (open-loop unstable) and leave
// the model's domain.
void BM_F16Rk4Step(benchmark::State& state) {
  const auto& op = operating();
  const double h = op.scenario.run->dt;
  auto f = [&](double tt, const Vec<17>& xx) { return op.scenario.model.deriv(tt, xx, op.trim.u); };
  StepStats stats;
  for (auto _ : state) {
    benchmark::DoNotOptimize(step(Integrator::RK4, f, 0.0, op.trim.x, h, true, stats));
  }
}
BENCHMARK(BM_F16Rk4Step);

// 60 s of C1 at dt = 0.01 s, nothing written: output_every is past the end, so only the first and last states
// are kept. The run includes its own trim solve (BM_F16Trim), so the real-time factor below is conservative.
void BM_F16ClosedLoop60s(benchmark::State& state) {
  sim::AircraftScenario s = operating().scenario;
  s.run->t_final = 60.0;
  s.run->output_every = 1 << 20;
  for (auto _ : state) {
    auto out = sim::run_aircraft(s);
    benchmark::DoNotOptimize(out);
  }
}
BENCHMARK(BM_F16ClosedLoop60s)->Unit(benchmark::kMillisecond);

// Trim and linearization at the S&L example condition (502 ft/s, sea level, ideal surfaces; linearize
// requires them).
const sim::AircraftScenario& ideal() {
  static const sim::AircraftScenario s =
      sim::load_aircraft_scenario(root / "scenarios/f16_trim_502_xcg035.toml");
  return s;
}

void BM_F16Trim(benchmark::State& state) {
  const auto& s = ideal();
  for (auto _ : state) {
    benchmark::DoNotOptimize(trim_level(s.model, s.trim));
  }
}
BENCHMARK(BM_F16Trim)->Unit(benchmark::kMicrosecond);

void BM_F16Linearize(benchmark::State& state) {
  const auto& s = ideal();
  const TrimResult tr = trim_level(s.model, s.trim);
  for (auto _ : state) {
    benchmark::DoNotOptimize(linearize(s.model, tr));
  }
}
BENCHMARK(BM_F16Linearize)->Unit(benchmark::kMicrosecond);

}  // namespace
