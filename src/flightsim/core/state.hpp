#pragma once

#include "flightsim/math/quat.hpp"

namespace flightsim::rb {

// Rigid-body state layout (PRD 4.4): [ p_N (3) | v_N (3) | q_NB (4) | omega_B (3) ].
inline constexpr int kPos = 0, kVel = 3, kQuat = 6, kOmega = 10, kN = 13;

// Segment accessors for any state with the rigid-body block first (N >= 13, e.g. the 17-state F-16).
#define FLIGHTSIM_RB_ACCESSOR(name, offset, len) \
  template <int N>                               \
  auto name(Vec<N>& x) {                         \
    static_assert(N >= kN);                      \
    return x.template segment<len>(offset);      \
  }                                              \
  template <int N>                               \
  auto name(const Vec<N>& x) {                   \
    static_assert(N >= kN);                      \
    return x.template segment<len>(offset);      \
  }
FLIGHTSIM_RB_ACCESSOR(pos, kPos, 3)
FLIGHTSIM_RB_ACCESSOR(vel, kVel, 3)
FLIGHTSIM_RB_ACCESSOR(quat, kQuat, 4)
FLIGHTSIM_RB_ACCESSOR(omega, kOmega, 3)
#undef FLIGHTSIM_RB_ACCESSOR

}  // namespace flightsim::rb
