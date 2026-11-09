#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/rigid_body.hpp"
#include "flightsim/core/state.hpp"
#include "flightsim/math/quat.hpp"

using namespace flightsim;

namespace {
Vec<13> state(const Vec3& p, const Vec3& v, const Vec4& q, const Vec3& w) {
  Vec<13> x;
  x << p, v, q, w;
  return x;
}
MassProps unit_body() { return MassProps{1.0, Mat3::Identity()}; }
}  // namespace

TEST(RigidBody, StateLayoutAndAccessors) {
  EXPECT_EQ(rb::kPos, 0);
  EXPECT_EQ(rb::kVel, 3);
  EXPECT_EQ(rb::kQuat, 6);
  EXPECT_EQ(rb::kOmega, 10);
  EXPECT_EQ(rb::kN, 13);
  Vec<14> x = Vec<14>::Zero();  // accessors accept any N >= 13
  for (int i = 0; i < 13; ++i) x(i) = i;
  EXPECT_EQ(rb::pos(x), Vec3(0, 1, 2));
  EXPECT_EQ(rb::vel(x), Vec3(3, 4, 5));
  EXPECT_EQ(rb::quat(x), Vec4(6, 7, 8, 9));
  EXPECT_EQ(rb::omega(x), Vec3(10, 11, 12));
  rb::omega(x) = Vec3(1, 2, 3);  // non-const accessors write through
  EXPECT_EQ(x.segment<3>(10), Vec3(1, 2, 3));
  const Vec<14>& cx = x;
  EXPECT_EQ(rb::vel(cx), Vec3(3, 4, 5));
}

TEST(RigidBody, GravityOnlyAcceleration) {
  const Vec3 g(0, 0, 9.80665);
  const Vec<13> d = rigid_body_deriv(state(Vec3::Zero(), Vec3(1, 2, 3), Vec4(1, 0, 0, 0), Vec3::Zero()),
                                     unit_body(), ForceMoment{}, g);
  EXPECT_EQ(d.segment<3>(rb::kVel), g);              // vdot = g_N exactly (m = 1, no forces)
  EXPECT_EQ(d.segment<3>(rb::kPos), Vec3(1, 2, 3));  // pdot = v
  EXPECT_EQ(d.segment<3>(rb::kOmega), Vec3::Zero());
}

TEST(RigidBody, BodyForceRotatedToNed) {
  const MassProps mp{2.0, Mat3::Identity()};
  ForceMoment fm;
  fm.F_B = Vec3(4.0, 0, 0);
  const Vec<13> x =
      state(Vec3::Zero(), Vec3::Zero(), quat_from_euler(0, 0, std::numbers::pi / 2), Vec3::Zero());
  const Vec<13> d = rigid_body_deriv(x, mp, fm, Vec3::Zero());
  // psi = 90 deg: body x -> East, so vdot_N = [0, F/m, 0]. 1e-15: half-angle trig roundoff (cos(pi/4)^2
  // etc.).
  EXPECT_LE((d.segment<3>(rb::kVel) - Vec3(0, 2.0, 0)).norm(), 1e-15);
  // NED-expressed force is not rotated.
  ForceMoment fn;
  fn.F_N = Vec3(4.0, 0, 0);
  EXPECT_LE((rigid_body_deriv(x, mp, fn, Vec3::Zero()).segment<3>(rb::kVel) - Vec3(2.0, 0, 0)).norm(), 0.0);
}

TEST(RigidBody, QuatKinematicsMatchesFormula) {
  const Vec4 q = quat_normalized(Vec4(0.9, 0.1, -0.3, 0.2));
  const Vec3 w(0.3, -0.2, 0.5);
  const Vec<13> d =
      rigid_body_deriv(state(Vec3::Zero(), Vec3::Zero(), q, w), unit_body(), ForceMoment{}, Vec3::Zero());
  // qdot = 1/2 q (x) [0, w_B]: body rates multiply on the RIGHT (PRD 4.2). Hand-expanded components, not
  // quat_mul:
  const Vec4 expected = 0.5 * Vec4(-q(1) * w.x() - q(2) * w.y() - q(3) * w.z(),  //
                                   q(0) * w.x() + q(2) * w.z() - q(3) * w.y(),   //
                                   q(0) * w.y() - q(1) * w.z() + q(3) * w.x(),   //
                                   q(0) * w.z() + q(1) * w.y() - q(2) * w.x());
  EXPECT_LE((d.segment<4>(rb::kQuat) - expected).norm(),
            1e-16);  // identical arithmetic up to ordering: a few 1e-17
  // The left-multiplied variant 1/2 [0,w] (x) q differs for a generic attitude, so this pin is not vacuous.
  EXPECT_GT((0.5 * quat_mul(Vec4(0, w.x(), w.y(), w.z()), q) - expected).norm(), 1e-2);
}

TEST(RigidBody, EulerEquationsGyroscopic) {
  const MassProps mp{1.0, Mat3(Eigen::Vector3d(1, 2, 3).asDiagonal())};
  const Vec<13> d = rigid_body_deriv(state(Vec3::Zero(), Vec3::Zero(), Vec4(1, 0, 0, 0), Vec3(1, 1, 1)), mp,
                                     ForceMoment{}, Vec3::Zero());
  // wdot = I^-1(-w x I w) with I = diag(1,2,3), w = (1,1,1): w x Iw = (1,-2,1) -> (-1, 1, -1/3). Hand
  // computed.
  EXPECT_LE((d.segment<3>(rb::kOmega) - Vec3(-1.0, 1.0, -1.0 / 3.0)).norm(), 1e-15);
}

TEST(RigidBody, RotorTermIncluded) {
  MassProps mp{1.0, Mat3(Eigen::Vector3d(1, 2, 3).asDiagonal())};
  mp.h_rotor = Vec3(160, 0, 0);
  const Vec3 w(0, 1, 0);
  const Vec<13> d = rigid_body_deriv(state(Vec3::Zero(), Vec3::Zero(), Vec4(1, 0, 0, 0), w), mp,
                                     ForceMoment{}, Vec3::Zero());
  // w x h = (0,0,-160), w x I w = 0 (parallel), so wdot_z = -(w x h)_z / Izz = 160/3. Hand computed.
  EXPECT_NEAR(d(rb::kOmega + 2), 160.0 / 3.0, 1e-13);
  EXPECT_NEAR(d(rb::kOmega + 0), 0.0, 1e-15);
}

TEST(RigidBody, AppliedMomentAndForceModels) {
  const Vec<13> x = state(Vec3(1, 2, 3), Vec3(4, 5, 6), Vec4(1, 0, 0, 0), Vec3::Zero());
  EXPECT_EQ(linear_drag(0.5)(0.0, x).F_N, Vec3(-2.0, -2.5, -3.0));
  EXPECT_EQ(quadratic_drag(0.1)(0.0, x).F_N, -0.1 * std::sqrt(77.0) * Vec3(4, 5, 6));
  EXPECT_EQ(spring(Vec3(1, 4, 9))(0.0, x).F_N, Vec3(-1, -8, -27));
  EXPECT_EQ(body_force(Vec3(1, 2, 3))(0.0, x).F_B, Vec3(1, 2, 3));
  EXPECT_EQ(body_torque(Vec3(1, 2, 3))(0.0, x).M_B, Vec3(1, 2, 3));
  ForceMoment a, b;
  a.F_N = Vec3(1, 0, 0);
  b.F_N = Vec3(0, 1, 0);
  b.M_B = Vec3(0, 0, 1);
  a += b;
  EXPECT_EQ(a.F_N, Vec3(1, 1, 0));
  EXPECT_EQ(a.M_B, Vec3(0, 0, 1));
}
