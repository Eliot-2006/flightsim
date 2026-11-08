#include <benchmark/benchmark.h>

#include <Eigen/Dense>

static void BM_Vec3Dot(benchmark::State& state) {
  const Eigen::Vector3d a(1.0, 2.0, 3.0);
  const Eigen::Vector3d b(4.0, 5.0, 6.0);
  for (auto _ : state) {
    double d = a.dot(b);
    benchmark::DoNotOptimize(d);
  }
}
BENCHMARK(BM_Vec3Dot);

BENCHMARK_MAIN();
