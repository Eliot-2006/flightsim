#pragma once

#include <functional>
#include <string>
#include <vector>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/rigid_body.hpp"

namespace flightsim::validation {

// A validation problem with a known (or fine-reference) solution.
struct CaseDef {
  std::string name;
  MassProps mp;
  Vec3 g_N;
  std::vector<ForceModel> forces;
  Vec<13> x0;
  double t_final;
  std::function<Vec<13>(double t)> exact;  // empty -> no full-state closed form
  std::function<double(const Vec<13>&)> energy =
      {};                                          // Optional case invariant; default is rotational energy.
  std::function<Vec3(double t)> exact_omega = {};  // R3: elliptic rates, attitude needs a fine reference
};

// Translational cases (PRD 8.2); parameters are fixed and documented in cases.cpp.
CaseDef case_T0();  // zero force (smoke)
CaseDef case_T1();  // constant acceleration (smoke, RK4 exact)
CaseDef case_T2();  // 3D linear drag + gravity
CaseDef case_T3();  // vertical quadratic drag from rest
CaseDef case_T4();  // anisotropic harmonic oscillator
CaseDef case_T5();  // rotating body-x force (catches C_NB vs C_BN)

CaseDef case_R1();
CaseDef case_R2();
CaseDef case_R3();
CaseDef case_R4();
CaseDef case_R5();

}  // namespace flightsim::validation
