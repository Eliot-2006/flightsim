#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <random>

#include "flightsim/math/quat.hpp"
#include "flightsim/math/units.hpp"

using namespace flightsim;

namespace {

Mat3 Rx(double a) { return Eigen::AngleAxisd(a, Vec3::UnitX()).toRotationMatrix(); }
Mat3 Ry(double a) { return Eigen::AngleAxisd(a, Vec3::UnitY()).toRotationMatrix(); }
Mat3 Rz(double a) { return Eigen::AngleAxisd(a, Vec3::UnitZ()).toRotationMatrix(); }

constexpr double kDeg = units::deg;

}  // namespace

TEST(Quat, HamiltonProductMatchesHandComputed) {
  const Vec4 i(0, 1, 0, 0), j(0, 0, 1, 0), k(0, 0, 0, 1);
  // Exact small-integer arithmetic: no rounding, so equality is exact.
  EXPECT_EQ(quat_mul(i, j), k);
  EXPECT_EQ(quat_mul(j, i), -k);
  EXPECT_EQ(quat_mul(i, i), Vec4(-1, 0, 0, 0));
  // (1+2i+3j+4k)(5+6i+7j+8k) = -60+12i+30j+24k (hand computed; integers, exact).
  EXPECT_EQ(quat_mul(Vec4(1, 2, 3, 4), Vec4(5, 6, 7, 8)), Vec4(-60, 12, 30, 24));
}

TEST(Quat, ConjugateInverse) {
  const Vec4 q = quat_normalized(Vec4(0.3, -0.5, 0.7, 0.2));
  const Vec4 p = quat_mul(q, quat_conj(q));
  // Tolerance: a handful of double roundoffs on unit-magnitude products (plan value).
  EXPECT_LE((p - Vec4(1, 0, 0, 0)).norm(), 1e-15);
}

TEST(Quat, RotateBasis90) {
  const double h = 90.0 * kDeg;
  // Tolerance 1e-15: cos(pi/2) = 6e-17 plus half-angle roundoff, a few ulp of unit values.
  EXPECT_LE((quat_rotate(quat_from_euler(0, 0, h), Vec3::UnitX()) - Vec3(0, 1, 0)).norm(),
            1e-15);  // psi=90: East
  EXPECT_LE((quat_rotate(quat_from_euler(0, h, 0), Vec3::UnitX()) - Vec3(0, 0, -1)).norm(),
            1e-15);  // theta=90: Up
  EXPECT_LE((quat_rotate(quat_from_euler(h, 0, 0), Vec3::UnitY()) - Vec3(0, 0, 1)).norm(),
            1e-15);  // phi=90: y -> Down
}

TEST(Quat, RotateMatchesSandwichProduct) {
  const Vec4 q = quat_from_euler(0.4, -0.3, 1.1);
  const Vec3 v(1.0, -2.0, 0.5);
  const Vec4 s =
      quat_mul(quat_mul(q, Vec4(0, v.x(), v.y(), v.z())), quat_conj(q));  // PRD 4.2: q (x) [0,v] (x) q*
  // Tolerance: two quaternion products vs one DCM product, a few ulp at magnitude ~2.
  EXPECT_LE((s.tail<3>() - quat_rotate(q, v)).norm(), 1e-15);
  EXPECT_NEAR(s(0), 0.0, 1e-15);
}

TEST(Quat, EulerSequence321MatchesElementaryRotations) {
  const double phi = 0.3, theta = -0.5, psi = 2.0;
  // PRD 4.1: C_NB = Rz(psi) Ry(theta) Rx(phi).
  const Mat3 expected = Rz(psi) * Ry(theta) * Rx(phi);
  // Tolerance: product of three rotation matrices, a few ulp.
  EXPECT_LE((dcm_from_quat(quat_from_euler(phi, theta, psi)) - expected).norm(), 1e-15);
}

TEST(Quat, DcmOrthonormal) {
  std::mt19937_64 rng(12345);
  std::normal_distribution<double> n01;
  double worst_orth = 0.0, worst_det = 0.0;
  for (int i = 0; i < 1000; ++i) {
    const Vec4 q = quat_normalized(Vec4(n01(rng), n01(rng), n01(rng), n01(rng)));
    const Mat3 C = dcm_from_quat(q);
    worst_orth = std::max(worst_orth, (C.transpose() * C - Mat3::Identity()).norm());
    worst_det = std::max(worst_det, std::abs(C.determinant() - 1.0));
  }
  // Tolerance 1e-14: plan value, ~50 ulp for a 3x3 product of unit-magnitude entries.
  EXPECT_LE(worst_orth, 1e-14);
  EXPECT_LE(worst_det, 1e-14);
}

TEST(Quat, EulerRoundTrip) {
  double worst = 0.0;
  std::vector<double> thetas = {-89.9, 89.9};
  for (int t = -80; t <= 80; t += 10) thetas.push_back(t);
  for (int p = -170; p <= 170; p += 10) {
    for (int s = -170; s <= 170; s += 10) {
      for (double t : thetas) {
        const Vec3 e(p * kDeg, t * kDeg, s * kDeg);
        const Vec3 r = euler_from_quat(quat_from_euler(e.x(), e.y(), e.z()));
        worst = std::max(worst, (r - e).cwiseAbs().maxCoeff());
      }
    }
  }
  // Tolerance 1e-12: plan value; error amplification 1/cos(theta) = 573 at 89.9 deg times ~1e-16.
  EXPECT_LE(worst, 1e-12);
}

TEST(Quat, GimbalLockDefined) {
  for (double theta : {90.0 * kDeg, -90.0 * kDeg}) {
    const Vec4 q = quat_from_euler(0.3, theta, 0.7);
    const Vec3 e = euler_from_quat(q);
    EXPECT_EQ(e.x(), 0.0);             // convention: phi = 0 at theta = +-90 deg
    EXPECT_NEAR(e.y(), theta, 1e-12);  // theta recovered (asin conditioning at the pole, plan value)
    // Attitude is still reproduced: psi absorbs the roll. Tolerance 1e-12 (plan value).
    EXPECT_LE((dcm_from_quat(quat_from_euler(e.x(), e.y(), e.z())) - dcm_from_quat(q)).norm(), 1e-12);
  }
}

TEST(Quat, ErrorAngleSignInvariant) {
  const Vec4 q = quat_from_euler(0.2, -0.4, 1.3);
  EXPECT_EQ(attitude_error(q, -q), 0.0);  // exact: q* (x) (-q) has exactly-zero vector part

  const Vec3 axis = Vec3(1, 2, 3).normalized();
  const Vec4 q1 = quat_mul(q, quat_from_rotvec(axis * (1.0 * kDeg)));
  // Tolerance 1e-12 (plan value): quaternion product roundoff (~1e-16) amplified by atan2 near 0.5 deg.
  EXPECT_NEAR(attitude_error(q, q1), 1.0 * kDeg, 1e-12);
  EXPECT_NEAR(attitude_error(q, -q1), 1.0 * kDeg, 1e-12);

  const Vec4 q2 = quat_mul(q, quat_from_rotvec(axis * 1e-9));
  // Tolerance 1e-15 abs (plan value): 1e-9 rad is resolved by atan2(|vec|,|w|), not acos(w) (which gives
  // ~1e-8).
  EXPECT_NEAR(attitude_error(q, q2), 1e-9, 1e-15);
}

TEST(Quat, RotvecRoundTrip) {
  const Vec3 axis = Vec3(1, -2, 3).normalized();
  for (double a : {0.0, 1e-10, 1e-5, 1.0, 3.1}) {
    const Vec3 rv = axis * a;
    const Vec4 q = quat_from_rotvec(rv);
    EXPECT_NEAR(q.norm(), 1.0, 1e-15) << a;
    // Tolerance 1e-14: atan2/trig roundoff, plus the small-angle Taylor branches (truncation < 1e-20 for
    // |rv|<1e-4).
    EXPECT_LE((rotvec_from_quat(q) - rv).norm(), 1e-14) << a;
    EXPECT_LE((rotvec_from_quat(-q) - rv).norm(), 1e-14) << a;  // w >= 0 canonicalization
    EXPECT_LE((quat_from_rotvec(rotvec_from_quat(q)) - q).norm(), 1e-14) << a;
  }
  // Zero rotation is the identity quaternion exactly.
  EXPECT_EQ(quat_from_rotvec(Vec3::Zero()), Vec4(1, 0, 0, 0));
  EXPECT_EQ(rotvec_from_quat(Vec4(1, 0, 0, 0)), Vec3::Zero());
}

TEST(Quat, EigenStorageOrderGuard) {
  const Vec4 q = quat_normalized(Vec4(0.9, 0.1, -0.3, 0.2));
  // Eigen::Quaterniond(w,x,y,z) is scalar-first at construction but stores (x,y,z,w); we only compare
  // rotations.
  const Mat3 ref = Eigen::Quaterniond(q(0), q(1), q(2), q(3)).toRotationMatrix();
  EXPECT_LE((dcm_from_quat(q) - ref).norm(), 1e-15);  // a few ulp
}
