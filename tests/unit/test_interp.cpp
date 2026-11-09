#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>

#include "flightsim/math/interp.hpp"

using namespace flightsim;

namespace {

Table1D t1() { return Table1D("t1", {0.0, 1.0, 2.0}, {0.0, 1.0, 4.0}); }

// v(i,j) = x_i^2 + y_j^2 sampled on x={0,1,2}, y={0,1,3}
Table2D t2() {
  Eigen::MatrixXd v(3, 3);
  v << 0, 1, 9,  //
      1, 2, 10,  //
      4, 5, 13;
  return Table2D("t2", {0.0, 1.0, 2.0}, {0.0, 1.0, 3.0}, v);
}

// Hand-computed values below are exact in binary FP, so tolerances are only rounding-level:
// 1e-12 absolute = ~4 orders above double eps for O(10) values.
constexpr double kTol = 1e-12;

}  // namespace

TEST(Interp, ExactAtBreakpoints) {
  const auto a = t1();
  EXPECT_EQ(a(0.0), 0.0);
  EXPECT_EQ(a(1.0), 1.0);
  EXPECT_EQ(a(2.0), 4.0);  // last breakpoint must be exact, not 1 ulp off
  const auto b = t2();
  const double xs[] = {0, 1, 2}, ys[] = {0, 1, 3};
  const double want[3][3] = {{0, 1, 9}, {1, 2, 10}, {4, 5, 13}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) EXPECT_EQ(b(xs[i], ys[j]), want[i][j]) << i << "," << j;
}

TEST(Interp, LinearBetweenBreakpoints) {
  EXPECT_NEAR(t1()(0.5), 0.5, kTol);
  EXPECT_NEAR(t1()(1.5), 2.5, kTol);
  // x=0.5,y=2 : lerp in y at x=0: 1+0.5*(9-1)=5 ; at x=1: 2+0.5*(10-2)=6 ; lerp in x -> 5.5
  EXPECT_NEAR(t2()(0.5, 2.0), 5.5, kTol);
}

TEST(Interp, ExtrapolatesFromBoundaryCell) {
  EXPECT_NEAR(t1()(3.0), 7.0, kTol);      // last cell slope 3: 4 + 3*1
  EXPECT_NEAR(t1()(-1.0), -1.0, kTol);    // first cell slope 1
  EXPECT_NEAR(t1()(102.0), 304.0, kTol);  // no distance limit: 4 + 3*100
  const auto b = t2();
  EXPECT_NEAR(b(1.0, 4.0), 14.0, kTol);   // y high: row x=1 slope (10-2)/2 = 4, 10 + 4
  EXPECT_NEAR(b(0.0, -1.0), -1.0, kTol);  // y low: row x=0 slope 1, 0 - 1
  EXPECT_NEAR(b(3.0, 4.0), 20.0, kTol);   // both axes high: cell x[1,2] y[1,3], bilinear s=2, t=1.5
  EXPECT_NEAR(b(-1.0, 0.0), -1.0, kTol);  // x low: column y=0 slope 1
}

TEST(Interp, MonitorFlagsOutOfDomain) {
  LookupMonitor mon;
  EXPECT_NEAR(t1()(1.5, &mon), 2.5, kTol);
  EXPECT_FALSE(mon.out_of_domain);
  EXPECT_TRUE(mon.first_violation.empty());

  t2()(0.5, 5.0, &mon);
  ASSERT_TRUE(mon.out_of_domain);
  EXPECT_NE(mon.first_violation.find("t2"), std::string::npos) << mon.first_violation;
  EXPECT_NE(mon.first_violation.find("y="), std::string::npos) << mon.first_violation;
  EXPECT_NE(mon.first_violation.find("[0,3]"), std::string::npos) << mon.first_violation;

  const std::string first = mon.first_violation;
  t1()(-4.0, &mon);  // a later violation must not overwrite the first
  EXPECT_EQ(mon.first_violation, first);

  mon.reset();
  EXPECT_FALSE(mon.out_of_domain);
  EXPECT_TRUE(mon.first_violation.empty());

  t1()(std::nan(""), &mon);  // NaN input is never "in domain"
  EXPECT_TRUE(mon.out_of_domain);

  EXPECT_NO_THROW(t1()(5.0, nullptr));  // null monitor allowed
}

TEST(Interp, MonitorRecordsCells) {
  LookupMonitor mon;
  t1()(0.5, &mon);  // record_cells off by default
  EXPECT_TRUE(mon.cells.empty());

  mon.record_cells = true;
  t1()(0.5, &mon);  // cell 0
  t1()(1.5, &mon);  // cell 1
  t1()(9.0, &mon);  // extrapolated: boundary cell 1
  ASSERT_EQ(mon.cells.size(), 3u);
  EXPECT_EQ(mon.cells[0][0], mon.cells[1][0]);  // same table id
  EXPECT_EQ(mon.cells[0][1], 0);
  EXPECT_EQ(mon.cells[1][1], 1);
  EXPECT_EQ(mon.cells[2][1], 1);

  mon.reset();
  EXPECT_TRUE(mon.cells.empty());
  EXPECT_TRUE(mon.record_cells);  // reset clears results, keeps configuration

  t2()(0.5, 2.0, &mon);  // cell (0,1) of a 2x2-cell grid: index i*(ny-1)+j = 1
  t2()(1.5, 0.5, &mon);  // cell (1,0): index 2
  ASSERT_EQ(mon.cells.size(), 2u);
  EXPECT_EQ(mon.cells[0][0], mon.cells[1][0]);
  EXPECT_EQ(mon.cells[0][1], 1);
  EXPECT_EQ(mon.cells[1][1], 2);
}

TEST(Interp, RejectsNonIncreasingBreakpoints) {
  EXPECT_THROW(Table1D("bad", {0.0, 1.0, 1.0}, {0.0, 1.0, 2.0}), std::invalid_argument);
  EXPECT_THROW(Table1D("bad", {0.0, 2.0, 1.0}, {0.0, 1.0, 2.0}), std::invalid_argument);
  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(2, 2);
  EXPECT_THROW(Table2D("bad", {0.0, 0.0}, {0.0, 1.0}, v), std::invalid_argument);
  EXPECT_THROW(Table2D("bad", {0.0, 1.0}, {1.0, 0.0}, v), std::invalid_argument);
}

TEST(Interp, RejectsShapeMismatch) {
  EXPECT_THROW(Table1D("bad", {0.0, 1.0}, {0.0, 1.0, 2.0}), std::invalid_argument);
  EXPECT_THROW(Table1D("bad", {0.0}, {0.0}), std::invalid_argument);  // need >= 1 cell
  EXPECT_THROW(Table2D("bad", {0.0, 1.0}, {0.0, 1.0, 2.0}, Eigen::MatrixXd::Zero(2, 2)),
               std::invalid_argument);
  EXPECT_THROW(Table2D("bad", {0.0, 1.0}, {0.0, 1.0}, Eigen::MatrixXd::Zero(3, 2)), std::invalid_argument);
}
