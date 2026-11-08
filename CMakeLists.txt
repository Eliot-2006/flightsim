cmake_minimum_required(VERSION 3.25)
project(flightsim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(Warnings)
include(Dependencies)

# Git sha injected at configure time (re-run cmake to refresh).
execute_process(COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE FLIGHTSIM_GIT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT FLIGHTSIM_GIT_SHA)
  set(FLIGHTSIM_GIT_SHA "unknown")
endif()

add_library(flightsim_core STATIC src/flightsim/version.cpp)
target_include_directories(flightsim_core PUBLIC src)
target_link_libraries(flightsim_core PUBLIC Eigen3::Eigen tomlplusplus::tomlplusplus)
target_compile_definitions(flightsim_core PRIVATE
  FLIGHTSIM_VERSION="${PROJECT_VERSION}" FLIGHTSIM_GIT_SHA="${FLIGHTSIM_GIT_SHA}")
flightsim_target_options(flightsim_core)

add_executable(flightsim apps/flightsim.cpp)
target_link_libraries(flightsim PRIVATE flightsim_core)
flightsim_target_options(flightsim)

enable_testing()
include(GoogleTest)

add_executable(unit_tests tests/unit/test_smoke.cpp)
target_link_libraries(unit_tests PRIVATE flightsim_core GTest::gtest_main)
flightsim_target_options(unit_tests)
gtest_discover_tests(unit_tests PROPERTIES LABELS unit)

add_executable(validation_tests tests/validation/test_smoke_validation.cpp)
target_link_libraries(validation_tests PRIVATE flightsim_core GTest::gtest_main)
flightsim_target_options(validation_tests)
gtest_discover_tests(validation_tests PROPERTIES LABELS validation)

add_executable(bench_flightsim bench/bench_smoke.cpp)
target_link_libraries(bench_flightsim PRIVATE flightsim_core benchmark::benchmark)
flightsim_target_options(bench_flightsim)
