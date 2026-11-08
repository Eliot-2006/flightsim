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
  std::function<Vec<13>(double t)> exact;  // empty -> no closed form
};

// Translational cases (PRD 8.2); parameters are fixed and documented in cases.cpp.
CaseDef case_T0();  // zero force (smoke)
CaseDef case_T1();  // constant acceleration (smoke, RK4 exact)
CaseDef case_T2();  // 3D linear drag + gravity
CaseDef case_T3();  // vertical quadratic drag from rest
CaseDef case_T4();  // anisotropic harmonic oscillator
CaseDef case_T5();  // rotating body-x force (catches C_NB vs C_BN)

}  // namespace flightsim::validation
