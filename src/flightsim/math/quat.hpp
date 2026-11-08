#pragma once

#include <Eigen/Dense>

namespace flightsim {

using Vec3 = Eigen::Vector3d;
using Vec4 = Eigen::Vector4d;
using Mat3 = Eigen::Matrix3d;
template <int N>
using Vec = Eigen::Matrix<double, N, 1>;

// Quaternions are Vec4 = [q0, q1, q2, q3], scalar first, Hamilton product (PRD 4.2).
// q_NB maps body components to NED components. Eigen::Quaterniond is never stored.
Vec4 quat_mul(const Vec4& a, const Vec4& b);
Vec4 quat_conj(const Vec4& q);
Vec4 quat_normalized(const Vec4& q);

Mat3 dcm_from_quat(const Vec4& q_NB);                 // C_NB (body -> NED); q should be unit
Vec3 quat_rotate(const Vec4& q_NB, const Vec3& v_B);  // C_NB v_B

// 3-2-1 Euler angles: C_NB = Rz(psi) Ry(theta) Rx(phi).
Vec4 quat_from_euler(double phi, double theta, double psi);
// Returns (phi, theta, psi). At theta = +-pi/2 (|C21|,|C22| < 1e-12) roll is undefined: phi = 0 and psi
// carries the whole heading/roll combination, so C_NB is still reproduced.
Vec3 euler_from_quat(const Vec4& q_NB);

// Rotation by |rv| about rv/|rv|. Taylor branch for |rv| < 1e-4 avoids 0/0.
Vec4 quat_from_rotvec(const Vec3& rv);
// Inverse of quat_from_rotvec; q and -q give the same result (shortest rotation, w >= 0).
Vec3 rotvec_from_quat(const Vec4& q);

// Sign-invariant rotation angle between two attitudes: 2*atan2(|vec(q_ref^* (x) q)|, |scalar|).
double attitude_error(const Vec4& q_ref, const Vec4& q);

Mat3 skew(const Vec3& v);

}  // namespace flightsim
