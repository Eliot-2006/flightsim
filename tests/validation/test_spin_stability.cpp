#include <gtest/gtest.h>

#include <Eigen/Eigenvalues>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <sstream>

#include "flightsim/analysis/jacobian.hpp"
#include "flightsim/validation/analytic.hpp"
#include "flightsim/validation/convergence.hpp"

using namespace flightsim;
using namespace flightsim::validation;

TEST(SpinStability, L1SpinStability) {
  const MassProps mp{1, Vec3(1, 2, 3).asDiagonal()};
  auto f = [&](const Vec3& w) -> Vec3 {
    Vec<13> x = Vec<13>::Zero();
    rb::quat(x) = quat_from_euler(.3, -.4, .5);
    rb::omega(x) = w;
    return rb::omega(rigid_body_deriv(x, mp, {}, Vec3::Zero()));
  };
  const Vec3 I(1, 2, 3);
  for (int i = 0; i < 3; ++i) {
    Vec3 w = Vec3::Zero();
    w(i) = 1.5;
    const Mat3 J = jacobian_central<3, 3>(f, w, Vec3::Ones());
    Eigen::EigenSolver<Mat3> solver(J);
    ASSERT_EQ(solver.info(), Eigen::Success);
    const int j = (i + 1) % 3, k = (i + 2) % 3;
    const double lambda2 = -2.25 * (I(i) - I(j)) * (I(i) - I(k)) / (I(j) * I(k));
    int nonzero = 0;
    bool positive = false;
    for (const auto& e : solver.eigenvalues()) {
      if (std::abs(e) < 1e-10) continue;  // Exact zero spin-axis mode, subtraction roundoff <1e-10.
      ++nonzero;
      // Task 2.3 analytic Euler spin eigenvalue formula, 1e-6 relative bound.
      EXPECT_LE(std::abs(e * e - lambda2) / std::abs(lambda2), 1e-6);
      positive = positive || e.real() > 1e-6;
    }
    EXPECT_EQ(nonzero, 2);
    EXPECT_EQ(positive, i == 1);
  }
}

TEST(Convergence, FineReferenceIndependentOfFirstOrderMethod) {
  auto c = case_T2();
  const auto truth = c.exact(c.t_final);
  c.exact = {};
  const auto rows = convergence_study(c, Integrator::Euler, .1, 8);
  for (const auto& row : rows) {
    const auto expected = state_error(integrate_case(c, Integrator::Euler, row.h).x_final, truth);
    // Fine RK4 error is near 1e-12; 1e-9 absolute is much smaller than same-Euler reference bias (~1e-3).
    EXPECT_NEAR(row.err.pos, expected.pos, 1e-9);
    EXPECT_NEAR(row.err.vel, expected.vel, 1e-9);
  }
}

TEST(Convergence, FineReferenceFourthOrderMethod) {
  auto c = case_T4();
  const auto truth = c.exact(c.t_final);
  c.exact = {};
  const auto rows = convergence_study(c, Integrator::RK4, .4, 8);
  const auto expected = state_error(integrate_case(c, Integrator::RK4, rows.back().h).x_final, truth);
  // PRD 6.3: fine-reference contamination <1e-3 of the finest method error.
  EXPECT_LE(std::abs(rows.back().err.pos - expected.pos), 1e-3 * expected.pos);
  EXPECT_LE(std::abs(rows.back().err.vel - expected.vel), 1e-3 * expected.vel);
}

TEST(Convergence, FineReferenceRefinesWhenInitialEstimateIsTooLarge) {
  auto c = case_T0();
  c.t_final = 1;
  c.exact = {};
  rb::pos(c.x0).setZero();
  rb::vel(c.x0).setZero();
  long evaluations = 0;
  c.forces = {[&evaluations](double t, const Vec<13>&) {
    ++evaluations;
    ForceMoment fm;
    fm.F_N.x() = std::cos(256 * std::numbers::pi * t);
    return fm;
  }};
  const auto rows = convergence_study(c, Integrator::RK4, 1, 4);
  // Initial /8 and /16 references alias this forcing; refinement resolves its zero final impulse.
  // Analytic method endpoint is (p,v)=(.5,1), exact endpoint is (0,0); trig roundoff bound 1e-10.
  EXPECT_NEAR(rows.back().err.pos, .5, 1e-10);
  EXPECT_NEAR(rows.back().err.vel, 1, 1e-10);
  EXPECT_GT(evaluations, 4 * (1 + 2 + 4 + 8 + 64 + 128));
  EXPECT_EQ(rows.back().nfev, 32);  // Work counts exclude reference construction.
}

TEST(Convergence, RejectsInvalidGridAndNonfiniteIntegration) {
  auto c = case_T0();
  EXPECT_THROW(convergence_study(c, Integrator::RK4, .1, 0), std::invalid_argument);
  EXPECT_THROW(integrate_case(c, Integrator::RK4, 0), std::invalid_argument);
  EXPECT_THROW(integrate_case(c, Integrator::RK4, std::nan("")), std::invalid_argument);
  c.forces = {[](double, const Vec<13>&) {
    ForceMoment f;
    f.F_N.setConstant(std::nan(""));
    return f;
  }};
  EXPECT_THROW(integrate_case(c, Integrator::RK4, .1), std::runtime_error);
}

TEST(Convergence, T4OrdersForApplicableIntegrators) {
  const auto c = case_T4();
  for (auto method : {Integrator::Euler, Integrator::RK4, Integrator::SymplecticEuler}) {
    const auto rows = convergence_study(c, method, method == Integrator::RK4 ? .025 : .00625, 8);
    std::vector<double> h, p, v;
    for (const auto& row : rows) {
      h.push_back(row.h);
      p.push_back(row.err.pos);
      v.push_back(row.err.vel);
      std::printf("[T4] %s h=%.12g pos=%.12g vel=%.12g\n", to_string(method), row.h, row.err.pos,
                  row.err.vel);
    }
    const bool fourth = method == Integrator::RK4;
    for (const auto& errors : {p, v}) {
      const auto fit = observed_order(h, errors, 1);
      // PRD 10.3, unchanged method order bounds.
      EXPECT_GE(fit.order, fourth ? 3.9 : .95) << to_string(method);
      EXPECT_LE(fit.order, fourth ? 4.1 : 1.05) << to_string(method);
    }
  }
}

TEST(Convergence, WritesStudyArtifactsAndRejectsInvalidInput) {
  const auto dir = std::filesystem::path(::testing::TempDir()) / "m2-study";
  std::filesystem::create_directories(dir);
  const auto input = dir / "study.toml";
  std::ofstream(input) << "[[studies]]\ncases=[\"T2\"]\nintegrators=[\"euler\"]\nh0=0.1\nlevels=8\n";
  EXPECT_NO_THROW(write_convergence_study(input, dir));
  std::ifstream csv(dir / "convergence.csv");
  ASSERT_TRUE(csv);
  std::string line;
  int lines = 0;
  while (std::getline(csv, line)) ++lines;
  EXPECT_EQ(lines, 9);  // Header + eight convergence levels.
  EXPECT_GT(std::filesystem::file_size(dir / "orders.json"), 0u);
  EXPECT_GT(std::filesystem::file_size(dir / "norm_conservation.csv"), 0u);
  for (const std::string text :
       {"[[studies]]\ncases=[\"T2\"]\nintegrators=[\"euler\"]\nh0=0.1\nlevels=8\nbogus=true\n",
        "[[studies]]\ncases=[\"nonsense\"]\nintegrators=[\"euler\"]\nh0=0.1\nlevels=8\n",
        "[[studies]]\ncases=[\"T2\"]\nintegrators=[\"euler\"]\nh0=0.1\nlevels=1\n",
        "[[studies]]\ncases=[\"T2\"]\nintegrators=[\"euler\"]\nh0=0.1\nlevels=8\nrenormalize=\"false\"\n"}) {
    std::ofstream(input) << text;
    EXPECT_THROW(write_convergence_study(input, dir), std::exception);
  }
}

TEST(Convergence, T4EnergyIncludesSpringPotential) {
  const auto c = case_T4();
  const auto rows = convergence_study(c, Integrator::Euler, .01, 2);
  const auto x = integrate_case(c, Integrator::Euler, .01).x_final;
  auto energy = [](const Vec<13>& y) {
    return .5 * (rb::vel(y).squaredNorm() + Vec3(1, 4, 9).dot(rb::pos(y).cwiseAbs2()));
  };
  // Same endpoint and arithmetic: 1e-14 relative invariant drift comparison (~45 ulp).
  EXPECT_NEAR(rows.front().err_E, std::abs(energy(x) / energy(c.x0) - 1), 1e-14);
}

TEST(Convergence, ConservationTracksTrajectoryMaxima) {
  auto c = case_R1();
  c.mp.inertia = Mat3::Identity();
  rb::omega(c.x0) = Vec3(0, 0, 1);
  c.t_final = 1;
  c.forces = {[](double t, const Vec<13>&) {
    ForceMoment fm;
    fm.M_B.z() = std::sin(2 * std::numbers::pi * t);
    return fm;
  }};
  const auto result = integrate_case(c, Integrator::RK4, .001);
  const double pi = std::numbers::pi;
  // w3=1+(1-cos(2*pi*t))/(2*pi), peaking at t=.5, returning to 1 at t=1.
  // RK4's sinusoidal quadrature bound O((2*pi)^4*h^4/180), propagated to E/H, is below 1e-10.
  EXPECT_NEAR(result.max_energy_drift, 2 / pi + 1 / (pi * pi), 1e-10);
  EXPECT_NEAR(result.max_momentum_drift, 1 / pi, 1e-10);
  EXPECT_LE(std::abs(rotational_energy(result.x_final, c.mp) / rotational_energy(c.x0, c.mp) - 1), 1e-10);
  EXPECT_GE(result.max_qnorm_drift, std::abs(rb::quat(result.x_final).norm() - 1));
}
