#include <cstring>
#include <iostream>

#include "flightsim/version.hpp"

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "--version") == 0) {
    std::cout << "flightsim " << flightsim::version_string() << '\n';
    return 0;
  }
  std::cerr << "usage: flightsim --version\n";
  return 2;
}
