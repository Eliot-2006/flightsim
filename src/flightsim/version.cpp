#include "flightsim/version.hpp"

#ifndef FLIGHTSIM_VERSION
#error "FLIGHTSIM_VERSION must be defined by the build"
#endif
#ifndef FLIGHTSIM_GIT_SHA
#error "FLIGHTSIM_GIT_SHA must be defined by the build"
#endif

namespace flightsim {

std::string version_string() { return std::string(FLIGHTSIM_VERSION) + " (" + FLIGHTSIM_GIT_SHA + ")"; }

}  // namespace flightsim
