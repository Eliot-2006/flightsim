#pragma once

#include <filesystem>

#include "flightsim/sim/runner.hpp"

namespace flightsim::sim {

// One row per sample. Columns carry units: t_s, pN/pE/pD_m, vN/vE/vD_mps, q0..q3, p/q/r_radps, then derived
// phi/theta/psi_rad (3-2-1 Euler), h_m (= -pD), V_mps/alpha_rad/beta_rad (air-relative = inertial; no wind in
// the generic body), E_J (kinetic + rotational + m g h), HN_{x,y,z}_kgm2ps (C_NB (I w + h_rotor)), qnorm.
// Numbers are printed with %.17g (round-trip exact), so identical runs give identical bytes.
void write_csv(const std::filesystem::path&, const RigidBodyScenario&, const RunOutput&);

// JSON: version, git sha, build type, compiler, scenario path + hash, resolved run parameters, stop reason,
// and the scenario file contents (all force parameters).
void write_metadata(const std::filesystem::path&, const RigidBodyScenario&, const RunOutput&);

}  // namespace flightsim::sim
