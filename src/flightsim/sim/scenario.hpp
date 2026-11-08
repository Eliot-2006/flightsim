#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "flightsim/core/forces.hpp"
#include "flightsim/core/integrators.hpp"
#include "flightsim/core/rigid_body.hpp"

namespace flightsim::sim {

struct RunConfig {
  Integrator integ;
  double dt, t_final;
  int output_every = 1;
  bool renormalize = true;
  bool stop_at_ground = false;  // stop at the first step boundary with h = -p_D <= 0
};

struct RigidBodyScenario {
  RunConfig run;
  MassProps mp;
  Vec3 g_N;
  std::vector<ForceModel> forces;
  Vec<13> x0;
  std::string source_path, source_hash;  // hash = "fnv1a64:<16 hex>" of the file bytes
};

// Schema = PRD 9.4: [sim], [body], [initial], optional [environment], optional [[forces]].
//   [sim]      integrator, dt_s, t_final_s, output_every (1), renormalize (true), stop_at_ground (false)
//   [body]     mass_kg, inertia_kgm2 (3x3 symmetric positive definite), rotor_angular_momentum_kgm2ps
//   ([0,0,0]) [initial]  position_m, velocity_mps, euler_deg (phi, theta, psi), omega_radps (body rates)
//   [environment]  gravity_mps2 (9.80665), applied along +z (down)
//   [[forces]] type = linear_drag{c_kgps} | quadratic_drag{k_kgpm} | spring{k_Npm[3]} |
//   body_force{force_N[3]} |
//              body_torque{torque_Nm[3]}
// Throws std::runtime_error naming the offending key path for: unknown keys, missing/mistyped keys, a dt that
// does not divide t_final (or the output interval), and non-positive-definite inertia.
RigidBodyScenario load_rigid_body_scenario(const std::filesystem::path&);

}  // namespace flightsim::sim
