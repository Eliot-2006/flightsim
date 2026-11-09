#include <gtest/gtest.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <toml++/toml.hpp>

#include "flightsim/analysis/flight_vars.hpp"
#include "flightsim/analysis/linearize.hpp"
#include "flightsim/analysis/modes.hpp"
#include "flightsim/math/units.hpp"
#include "flightsim/sim/aircraft_scenario.hpp"

using namespace flightsim;
namespace {
const auto root = std::filesystem::path(FLIGHTSIM_SOURCE_DIR);
auto scenario() { return sim::load_aircraft_scenario(root / "scenarios/f16_trim_502_xcg030.toml"); }
auto reference() { return toml::parse_file((root / "data/f16/reference_values.toml").string()); }
void matrix_check(const Eigen::MatrixXd& ours, const toml::array& ref, double relative) {
  for (Eigen::Index i = 0; i < ours.rows(); ++i) {
    for (Eigen::Index j = 0; j < ours.cols(); ++j) {
      const double p =
          ref[static_cast<size_t>(i)].as_array()->get(static_cast<size_t>(j))->value<double>().value();
      // Per-entry last printed digit from LN1 arrays, or Task 5.2 relative bound.
      // TOML preserves literals in source; source region gives the decimal precision.
      const auto node = ref[static_cast<size_t>(i)].as_array()->get(static_cast<size_t>(j));
      const auto region = node->source();
      std::ifstream file(root / "data/f16/reference_values.toml");
      std::string line;
      for (uint32_t n = 0; n < region.begin.line; ++n) {
        std::getline(file, line);
      }
      const std::string literal =
          line.substr(region.begin.column - 1, region.end.column - region.begin.column);
      const auto dot = literal.find('.');
      const double digit = dot == std::string::npos
                               ? (p == 0 ? 0 : 1)
                               : std::pow(10.0, -static_cast<double>(literal.size() - dot - 1));
      const double difference = ours(i, j) - p;
      std::cout << "LN1 matrix [" << i << ',' << j << "] ours=" << ours(i, j) << " published=" << p
                << " difference=" << difference << '\n';
      // Exact structural zeros tolerate central-difference roundoff, observed <1e-10, 10x safety.
      EXPECT_LE(std::abs(difference), std::max({relative * std::abs(p), digit, 1e-9}));
    }
  }
}
}  // namespace
TEST(F16Linear, LN1Level30Eigenvalues) {
  const auto s = scenario();
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  const auto ref = reference();
  // Published model has NO power state: identify_modes uses [V,alpha,theta,q] with P held at trim.
  const auto modes = identify_modes(lm);
  ASSERT_EQ(modes.size(), 5u);
  for (const auto& mode : modes) {
    const std::string key = mode.name == "roll" ? "roll_subsidence" : mode.name;
    auto node = ref["LN1"]["level_30"][key]["value"];
    std::complex<double> p;
    if (const auto* a = node.as_array()) {
      p = {(*a)[0].value<double>().value(), (*a)[1].value<double>().value()};
    } else {
      p = {node.value<double>().value(), 0};
    }
    std::cout << std::setprecision(12) << "LN1 " << key << " ours=" << mode.eig << " published=" << p
              << " difference=" << mode.eig - p << '\n';
    // PRD 8.5 LN1 acceptance: eigenvalues within 1 % (Task 5.2: or 1 in the last printed digit), per
    // real/imaginary component. This test enforces the tighter regression bound max(2 x observed, 1 in the
    // last printed digit) as max(2.5e-4 relative, digit): the only component outside one printed digit is
    // the short-period imaginary part (1.72e-4 against 1e-4, i.e. 1.16e-4 relative); all others are within
    // their digit (phugoid real 2.95e-5 vs 1e-4, roll 8.26e-4 vs 1e-3, dutch-roll imaginary 3.34e-4 vs 1e-2).
    constexpr double kEigRelative = 2.5e-4;
    EXPECT_LE(std::abs(mode.eig.real() - p.real()),
              std::max(kEigRelative * std::abs(p.real()), key == "roll_subsidence" ? .001 : .0001));
    EXPECT_LE(std::abs(mode.eig.imag() - p.imag()),
              std::max(kEigRelative * std::abs(p.imag()), key == "dutch_roll" ? .01 : .0001));
    EXPECT_DOUBLE_EQ(mode.wn, std::abs(mode.eig));
    EXPECT_DOUBLE_EQ(mode.zeta, -mode.eig.real() / mode.wn);
    EXPECT_DOUBLE_EQ(mode.tau, -1 / mode.eig.real());
  }
}
TEST(F16Linear, LN1LateralMatrix) {
  auto s = scenario();
  auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  const auto ref = reference();
  Eigen::Matrix4d a;
  const std::array<int, 4> ix{2, 3, 6, 8};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      a(i, j) = lm.A(ix[static_cast<std::size_t>(i)], ix[static_cast<std::size_t>(j)]);
    }
  }
  // PRD 8.5 / Task 5.2 acceptance: 1 % or one last printed digit. Enforced here: the regression bound
  // max(2 x observed, digit) as max(2e-3 relative, digit). Entries outside their printed digit reach at most
  // 9.23e-4 relative (r-row, p-column -0.026358, difference 2.43e-5); the others are at most 2.44e-4.
  matrix_check(a, *ref["LN1"]["level_30"]["lateral_textbook"]["A"]["value"].as_array(), 2e-3);
}
TEST(F16Linear, LN1LongitudinalTutorialMatrix) {
  auto s = scenario();
  auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  const auto ref = reference();
  Eigen::Matrix4d a;
  Eigen::Matrix<double, 4, 2> b;
  const std::array<int, 4> ix{0, 1, 4, 7};
  for (int i = 0; i < 4; ++i) {
    const double row = i == 0 ? 1 / units::ft : 1;
    for (int j = 0; j < 4; ++j) {
      a(i, j) = row * lm.A(ix[static_cast<std::size_t>(i)], ix[static_cast<std::size_t>(j)]) *
                (j == 0 ? units::ft : 1);
    }
    // Tutorial-computed B uses tabulated engine: eliminate steady Pdot=0 for throttle, not free power.
    b(i, 0) = row * (lm.B(ix[static_cast<std::size_t>(i)], 0) -
                     lm.A(ix[static_cast<std::size_t>(i)], 12) * lm.B(12, 0) / lm.A(12, 12));
    b(i, 1) = row * lm.B(ix[static_cast<std::size_t>(i)], 1) * units::deg;
  }
  std::cout << "LN1 tutorial-computed (single-lineage secondary source)\n";
  // Task 5.2 acceptance for this rounded tutorial-computed matrix: 2 % or one last printed digit. The entries
  // are rounded to display precision (2-3 significant figures), and every observed difference is within one
  // printed digit (largest: B alpha-row -0.0021/-0.0022 at 0.50/0.49 digit, A[alpha,V] -0.0003 at 0.46
  // digit), so the regression bound max(2 x observed, digit) is the last printed digit alone.
  matrix_check(a, *ref["LN1"]["level_30"]["longitudinal_tutorial"]["A"]["value"].as_array(), 0);
  matrix_check(b, *ref["LN1"]["level_30"]["longitudinal_tutorial"]["B"]["value"].as_array(), 0);
}
TEST(F16Linear, LN2QuadraticScaling) {
  auto s = scenario();
  // Interior thrust-table altitude avoids the initial lift-loss transient below sea level.
  s.trim.h = 100 * units::ft;
  ASSERT_FALSE(s.model.actuators);
  ASSERT_TRUE(s.model.include_rotor);
  // Orchestrator LN2 ruling: the beta=0 fold is smooth only when every alpha-row intercept is zero.
  for (const auto* table : {&s.model.data.cl, &s.model.data.cn}) {
    ASSERT_EQ(table->y().front(), 0);
    for (double alpha : table->x()) {
      ASSERT_EQ((*table)(alpha, 0), 0) << table->name() << " alpha=" << alpha;
    }
  }
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  // Ideal surfaces and constant input hold da=dr=0 exactly throughout both runs.
  ASSERT_EQ(tr.u.da, 0);
  ASSERT_EQ(tr.u.dr, 0);
  const auto big = elevator_response(s.model, tr, lm, -0.5 * units::deg);
  const auto small = elevator_response(s.model, tr, lm, -0.05 * units::deg);
  for (const auto* response : {&big, &small}) {
    // The zero-multiplier exclusion is valid only with exactly zero effective lateral inputs at EVERY stage.
    EXPECT_EQ(response->max_abs_aileron, 0.0);
    EXPECT_EQ(response->max_abs_rudder, 0.0);
    std::set<std::string> changed;
    for (const auto& event : response->cell_changes) {
      std::cout << "LN2 " << event.table << " cells " << event.from << "->" << event.to << " at t=" << event.t
                << " ignored=" << event.ignored << ": " << event.reason << '\n';
      EXPECT_TRUE(event.ignored);
      changed.insert(event.table);
      if (event.table == s.model.data.cl.name() || event.table == s.model.data.cn.name()) {
        const auto& table = event.table == s.model.data.cl.name() ? s.model.data.cl : s.model.data.cn;
        EXPECT_EQ(event.from, event.to);  // same alpha row and first |beta| cell across the sign fold
        EXPECT_EQ(event.to % static_cast<int>(table.y().size() - 1), 0);
      } else {
        EXPECT_TRUE(event.table == s.model.data.dlda.name() || event.table == s.model.data.dldr.name() ||
                    event.table == s.model.data.dnda.name() || event.table == s.model.data.dndr.name());
      }
    }
    const std::set<std::string> expected{s.model.data.cl.name(),   s.model.data.cn.name(),
                                         s.model.data.dlda.name(), s.model.data.dldr.name(),
                                         s.model.data.dnda.name(), s.model.data.dndr.name()};
    EXPECT_EQ(changed, expected);
  }
  EXPECT_FALSE(big.out_of_domain) << big.first_violation;
  EXPECT_FALSE(small.out_of_domain) << small.first_violation;
  EXPECT_FALSE(big.cells_changed);
  EXPECT_FALSE(small.cells_changed);
  const double ratio = big.max_error / small.max_error;
  std::cout << std::setprecision(12) << "LN2 errors=" << big.max_error << ',' << small.max_error
            << " ratio=" << ratio << '\n';
  // PRD LN2: tenfold input gives quadratic error within [70,130].
  EXPECT_GE(ratio, 70);
  EXPECT_LE(ratio, 130);
}
TEST(F16Linear, LN3StepSweep) {
  auto s = scenario();
  auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  for (const auto& p : linearization_sweep(s.model, tr)) {
    std::cout << "LN3 eps=" << p.eps << " max_relative_entry_change=" << p.change << '\n';
  }
  const auto a = linearize(s.model, tr, 6e-6);
  for (double eps : {6e-5, 6e-7}) {
    auto b = linearize(s.model, tr, eps);
    const double change = std::max((a.A - b.A).cwiseAbs().maxCoeff() / a.A.cwiseAbs().maxCoeff(),
                                   (a.B - b.B).cwiseAbs().maxCoeff() / a.B.cwiseAbs().maxCoeff());
    std::cout << "LN3 chosen=6e-6 neighbor=" << eps << " change=" << change << '\n';
    // Task 5.2 flat-region criterion, unchanged.
    EXPECT_LE(change, 1e-5);
  }
}

TEST(F16Linear, RejectsFailedTrimAndBadEpsilon) {
  const auto s = scenario();
  EXPECT_THROW(linearize(s.model, TrimResult{}), std::invalid_argument);
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  EXPECT_THROW(linearize(s.model, tr, 0), std::invalid_argument);
}

TEST(F16Linear, ClimbingResponseRemovesNominalAltitudeTravel) {
  auto s = scenario();
  s.trim.gamma = 5 * units::deg;
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  const auto response = elevator_response(s.model, tr, lm, -.5 * units::deg);
  ASSERT_GT(response.points.size(), 1u);
  const auto& point = response.points[1];
  // First .01s: altitude perturbation is <1e-5m from the elevator, versus .133m nominal climb.
  // Bound is the observed 8.1e-6m response rounded up, not a trajectory error tolerance.
  EXPECT_LT(std::abs(point.nonlinear(11)), 1e-5);
}

TEST(F16Linear, LN2RejectsActiveCellChanges) {
  auto s = scenario();
  s.trim.h = 100 * units::ft;
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  const auto response = elevator_response(s.model, tr, lm, -5 * units::deg);
  EXPECT_TRUE(response.cells_changed);
  // A larger step crosses active alpha cells. The exclusions must not mask these kinks.
  for (const auto* table : {&s.model.data.cx, &s.model.data.cm, &s.model.data.cl, &s.model.data.cn}) {
    bool found = false;
    for (const auto& event : response.cell_changes) {
      if (event.table == table->name() && event.from != event.to) {
        EXPECT_FALSE(event.ignored) << event.table;
        found = true;
      }
    }
    EXPECT_TRUE(found) << table->name();
  }
}

TEST(F16Linear, LN2EffectivenessChangesRequireZeroInputs) {
  auto s = scenario();
  s.trim.h = 100 * units::ft;
  auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  tr.u.da = .01 * units::deg;
  const auto response = elevator_response(s.model, tr, lm, -.5 * units::deg);
  EXPECT_EQ(response.max_abs_aileron, tr.u.da);
  EXPECT_TRUE(response.cells_changed);
  bool found = false;
  for (const auto& event : response.cell_changes) {
    if (event.table == s.model.data.dlda.name() || event.table == s.model.data.dldr.name() ||
        event.table == s.model.data.dnda.name() || event.table == s.model.data.dndr.name()) {
      EXPECT_FALSE(event.ignored) << event.table;
      found = true;
    }
  }
  EXPECT_TRUE(found);
}

TEST(F16Linear, LN2RejectsThrustCellChanges) {
  auto s = scenario();
  // Starting on an interior altitude breakpoint makes the initial lift transient cross a thrust cell.
  s.trim.h = 10000 * units::ft;
  const auto tr = trim_level(s.model, s.trim);
  ASSERT_TRUE(tr.ok);
  const auto response = elevator_response(s.model, tr, linearize(s.model, tr), -.5 * units::deg);
  EXPECT_TRUE(response.cells_changed);
  bool found = false;
  for (const auto& change : response.cell_changes) {
    if (change.table == s.model.data.thrust_mil.name()) {
      EXPECT_FALSE(change.ignored);
      found = true;
    }
  }
  EXPECT_TRUE(found);
}

TEST(F16Linear, LN2FoldExclusionRequiresEveryZeroIntercept) {
  const auto s = scenario();
  auto spec = s.trim;
  spec.h = 100 * units::ft;
  const auto tr = trim_level(s.model, spec);
  ASSERT_TRUE(tr.ok);
  const auto lm = linearize(s.model, tr);
  for (bool change_cl : {true, false}) {
    auto m = s.model;
    auto& table = change_cl ? m.data.cl : m.data.cn;
    Eigen::MatrixXd values(static_cast<Eigen::Index>(table.x().size()),
                           static_cast<Eigen::Index>(table.y().size()));
    for (Eigen::Index i = 0; i < values.rows(); ++i) {
      for (Eigen::Index j = 0; j < values.cols(); ++j) {
        values(i, j) = table(table.x()[static_cast<std::size_t>(i)], table.y()[static_cast<std::size_t>(j)]);
      }
    }
    // A nonzero intercept at even one remote alpha row invalidates the exact smooth-fold proof.
    values(0, 0) = 1e-12;
    table = Table2D(table.name(), table.x(), table.y(), values);
    const auto response = elevator_response(m, tr, lm, -.5 * units::deg);
    EXPECT_TRUE(response.cells_changed);
    bool found = false;
    for (const auto& event : response.cell_changes) {
      if (event.table == table.name() && event.from == event.to) {
        EXPECT_FALSE(event.ignored);
        found = true;
      }
    }
    EXPECT_TRUE(found);
  }
}
