#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <toml++/toml.hpp>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/analysis/trim.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;
namespace {
F16Model model(double cg) {
  F16Model m{load_f16_data(std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "data/f16/f16_sl.toml")};
  m.xcg = cg;
  m.grav.g0 = m.data.g_sl;
  return m;
}
void check(const char* name, double cg) {
  auto m = model(cg);
  LookupMonitor monitor;
  m.monitor = &monitor;
  const auto tr = trim_level(m, {502 * units::ft, 0, 0, cg});
  ASSERT_TRUE(tr.ok) << tr.message;
  monitor.reset();
  m.evaluate(0, tr.x, tr.u);
  EXPECT_FALSE(monitor.out_of_domain) << monitor.first_violation;
  const auto ref = toml::parse_file(std::string(FLIGHTSIM_SOURCE_DIR) + "/data/f16/reference_values.toml");
  const auto& t = *ref["TR1"][name].as_table();
  for (const auto& [key, value] : std::vector<std::pair<const char*, double>>{
           {"alpha", tr.alpha}, {"throttle", tr.throttle}, {"elevator", tr.de / units::deg}}) {
    const double published = t[key]["value"].value<double>().value();
    // PRD 8.5 TR1 acceptance: within 1 in the last printed digit or 0.5 % relative, whichever is larger.
    // This test enforces the tighter regression bound max(2 x observed, 1 in the last printed digit). Every
    // observed difference is below one printed digit (largest: level_30 elevator, 4.76e-4 deg against a
    // 1e-3 digit, i.e. 0.48 digit), so the bound is the last printed digit itself. The 0.5 % term would let a
    // 17x larger error through (worst observed relative difference 2.8e-4, level_38 throttle).
    const double digit = std::string(key) == "alpha"      ? 1e-5
                         : std::string(key) == "throttle" ? 1e-4
                         : cg == 0.30                     ? 1e-3
                                                          : 1e-4;
    const double difference = value - published;
    std::cout << std::setprecision(12) << "TR1 " << name << ' ' << key << " ours=" << value
              << " published=" << published << " difference=" << difference << '\n';
    EXPECT_LE(std::abs(difference), digit);
  }
}
}  // namespace
TEST(F16Trim, TR1Level35) { check("level_35", 0.35); }
TEST(F16Trim, TR1Level30) { check("level_30", 0.30); }
TEST(F16Trim, TR1Level38) { check("level_38", 0.38); }
TEST(F16Trim, TR2Stationarity) {
  auto m = model(0.35);
  const auto tr = trim_level(m, {502 * units::ft, 0, 0, 0.35});
  ASSERT_TRUE(tr.ok);
  LookupMonitor monitor;
  m.monitor = &monitor;
  Vec<17> x = tr.x;
  for (int i = 0; i < 1000; ++i) {
    x = rk4_step<17>([&](double t, const Vec<17>& y) { return m.deriv(t, y, tr.u); }, i * 0.01, x, 0.01);
    rb::quat(x) = quat_normalized(rb::quat(x));
  }
  const Vec<13> d = x13_from_x17(m, x) - x13_from_x17(m, tr.x);
  std::cout << "TR2 dV_fps=" << d(0) / units::ft << " dalpha=" << d(1) << " dtheta=" << d(4)
            << " dh_ft=" << d(11) / units::ft << '\n';
  // PRD TR2 stationarity bounds, unchanged.
  EXPECT_LT(std::abs(d(0)) / units::ft, 1e-3);
  EXPECT_LT(std::abs(d(1)), 1e-6);
  EXPECT_LT(std::abs(d(4)), 1e-6);
  EXPECT_LT(std::abs(d(11)) / units::ft, 1e-2);
  std::cout << "TR2 sea-level roundoff lookup: " << monitor.first_violation << '\n';
}
