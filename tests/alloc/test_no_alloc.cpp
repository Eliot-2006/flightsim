// PRD 3.1/11: the step loop performs no heap allocation. A global operator new counter is armed around the
// measured region; gtest's own allocations happen outside it.
#include <gtest/gtest.h>

#include <atomic>
#include <cstdlib>
#include <new>
#include <vector>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/math/quat.hpp"
#include "flightsim/math/units.hpp"

namespace {
std::atomic<bool> g_armed{false};
std::atomic<long> g_allocs{0};
}  // namespace

void* operator new(std::size_t n) {
  if (g_armed.load(std::memory_order_relaxed)) g_allocs.fetch_add(1, std::memory_order_relaxed);
  if (void* p = std::malloc(n ? n : 1)) return p;
  throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

namespace {
using namespace flightsim;

struct Counting {
  Counting() {
    g_allocs = 0;
    g_armed = true;
  }
  ~Counting() { g_armed = false; }
  long count() const { return g_allocs.load(); }
};

TEST(NoAlloc, CounterSeesHeapAllocation) {
  Counting c;
  auto* v = new std::vector<double>(100);
  const long seen = c.count();
  delete v;
  EXPECT_GE(seen, 1);  // proves the counter is live, so a zero below is meaningful
}

TEST(NoAlloc, Rk4StepLoopAllocatesNothing) {
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
  for (int i = 0; i < 100; ++i) x = step(Integrator::RK4, f, i * h, x, h, true, stats);  // warmup
  Counting c;
  for (int i = 0; i < 10000; ++i) x = step(Integrator::RK4, f, i * h, x, h, true, stats);
  EXPECT_EQ(c.count(), 0);
  EXPECT_TRUE(x.allFinite());
}
}  // namespace
