// PRD 11: rigid-body step cost per integrator, and a full R3 run per integrator.
#include <benchmark/benchmark.h>

#include <vector>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/validation/cases.hpp"
#include "flightsim/validation/convergence.hpp"

namespace {
using namespace flightsim;

constexpr Integrator kIntegrators[] = {Integrator::Euler, Integrator::RK4, Integrator::RKMK4};

// Gravity + quadratic drag, tumbling body: the PRD 11 "rigid-body RK4 step" workload.
void BM_Step(benchmark::State& state) {
  const Integrator integ = kIntegrators[state.range(0)];
  const MassProps mp{2.0, Vec3(1, 2, 3).asDiagonal()};
  const std::vector<ForceModel> forces{quadratic_drag(0.01)};
  const Vec3 g(0, 0, units::g0);
  auto f = [&](double t, const Vec<13>& x) {
    ForceMoment fm;
    for (const auto& force : forces) fm += force(t, x);
    return rigid_body_deriv(x, mp, fm, g);
  };
  Vec<13> x = Vec<13>::Zero();
  x.segment<4>(rb::kQuat) = Vec4(1, 0, 0, 0);
  x.segment<3>(rb::kVel) = Vec3(30, 5, -10);
  x.segment<3>(rb::kOmega) = Vec3(0.5, 0.1, 1.0);
  StepStats stats;
  const double h = 1e-3;
  double t = 0;
  for (auto _ : state) {
    x = step(integ, f, t, x, h, true, stats);
    t += h;
    benchmark::DoNotOptimize(x);
  }
  state.SetLabel(to_string(integ));
}
BENCHMARK(BM_Step)->DenseRange(0, 2);

// R3 (torque-free asymmetric body) for 10 s at h = 0.01 s: 1000 steps.
void BM_R3Run10s(benchmark::State& state) {
  const Integrator integ = kIntegrators[state.range(0)];
  validation::CaseDef c = validation::case_R3();
  c.t_final = 10.0;
  for (auto _ : state) {
    auto r = validation::integrate_case(c, integ, 0.01);
    benchmark::DoNotOptimize(r);
  }
  state.SetLabel(to_string(integ));
}
BENCHMARK(BM_R3Run10s)->DenseRange(0, 2)->Unit(benchmark::kMicrosecond);

}  // namespace

BENCHMARK_MAIN();
