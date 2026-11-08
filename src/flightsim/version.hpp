#pragma once

#include <string>

namespace flightsim {

// "<version> (<git sha>)"; the sha is injected at configure time.
std::string version_string();

}  // namespace flightsim
