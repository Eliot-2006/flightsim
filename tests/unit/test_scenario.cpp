#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "flightsim/sim/scenario.hpp"

using namespace flightsim;
using namespace flightsim::sim;

namespace {

const char* kBase = R"(
[sim]
integrator = "rk4"
dt_s = 0.01
t_final_s = 1.0
[body]
mass_kg = 2.0
inertia_kgm2 = [[1,0,0],[0,2,0],[0,0,3]]
[initial]
position_m = [0,0,0]
velocity_mps = [1,2,3]
euler_deg = [0,0,0]
omega_radps = [0,0,0]
)";

std::filesystem::path write_tmp(const std::string& name, const std::string& text) {
  const auto p = std::filesystem::path(::testing::TempDir()) / name;
  std::ofstream(p) << text;
  return p;
}

std::string replace(std::string s, const std::string& from, const std::string& to) {
  const auto i = s.find(from);
  EXPECT_NE(i, std::string::npos) << from;
  return s.replace(i, from.size(), to);
}

// Loads text and returns the exception message ("" if it loaded).
std::string load_error(const std::string& text) {
  try {
    load_rigid_body_scenario(write_tmp("scenario_test.toml", text));
  } catch (const std::runtime_error& e) {
    return e.what();
  }
  return "";
}

}  // namespace

TEST(Scenario, BaseScenarioLoads) {
  EXPECT_EQ(load_error(kBase), "");
  const auto s = load_rigid_body_scenario(write_tmp("base.toml", kBase));
  EXPECT_EQ(s.run.integ, Integrator::RK4);
  EXPECT_EQ(s.run.dt, 0.01);
  EXPECT_EQ(s.run.output_every, 1);
  EXPECT_TRUE(s.run.renormalize);
  EXPECT_FALSE(s.run.stop_at_ground);
  EXPECT_EQ(s.mp.mass, 2.0);
  EXPECT_EQ(s.mp.inertia(2, 2), 3.0);
  EXPECT_EQ(s.g_N, Vec3(0, 0, 9.80665));  // default gravity
  EXPECT_EQ(s.x0.segment<3>(3), Vec3(1, 2, 3));
  EXPECT_EQ(s.x0.segment<4>(6), Vec4(1, 0, 0, 0));
  EXPECT_EQ(s.source_hash.substr(0, 8), "fnv1a64:");
}

TEST(Scenario, RejectsDtNotDividingTFinal) {
  const std::string e = load_error(replace(kBase, "dt_s = 0.01", "dt_s = 0.3"));
  EXPECT_NE(e.find("sim.dt_s"), std::string::npos) << e;
  EXPECT_NE(e.find("divide"), std::string::npos) << e;
  // The output interval must divide the run too (final sample lands on t_final): 100 steps, every 3.
  const std::string e2 = load_error(replace(kBase, "dt_s = 0.01", "dt_s = 0.01\noutput_every = 3"));
  EXPECT_NE(e2.find("sim.output_every"), std::string::npos) << e2;
  EXPECT_EQ(load_error(replace(kBase, "dt_s = 0.01", "dt_s = 0.01\noutput_every = 5")), "");
  EXPECT_NE(load_error(replace(kBase, "dt_s = 0.01", "dt_s = -0.01")).find("sim.dt_s"), std::string::npos);
}

TEST(Scenario, RejectsUnknownKey) {
  EXPECT_NE(load_error(replace(kBase, "dt_s = 0.01", "dt_s = 0.01\nbogus = 1")).find("sim.bogus"),
            std::string::npos);
  EXPECT_NE(load_error(std::string(kBase) + "[nonsense]\nx = 1\n").find("nonsense"), std::string::npos);
  EXPECT_NE(load_error(std::string(kBase) + "[[forces]]\ntype = \"linear_drag\"\nc_kgps = 1.0\nextra = 2\n")
                .find("forces[0].extra"),
            std::string::npos);
}

TEST(Scenario, RejectsSingularInertia) {
  const std::string e = load_error(replace(kBase, "[[1,0,0],[0,2,0],[0,0,3]]", "[[1,0,0],[0,2,0],[0,0,0]]"));
  EXPECT_NE(e.find("body.inertia_kgm2"), std::string::npos) << e;
  // Not symmetric, and not positive definite (negative eigenvalue).
  EXPECT_NE(
      load_error(replace(kBase, "[[1,0,0],[0,2,0],[0,0,3]]", "[[1,1,0],[0,2,0],[0,0,3]]")).find("symmetric"),
      std::string::npos);
  EXPECT_NE(
      load_error(replace(kBase, "[[1,0,0],[0,2,0],[0,0,3]]", "[[1,0,0],[0,-2,0],[0,0,3]]")).find("positive"),
      std::string::npos);
}

TEST(Scenario, RejectsMissingField) {
  EXPECT_NE(load_error(replace(kBase, "mass_kg = 2.0\n", "")).find("body.mass_kg"), std::string::npos);
  EXPECT_NE(load_error(replace(kBase, "omega_radps = [0,0,0]", "")).find("initial.omega_radps"),
            std::string::npos);
  EXPECT_NE(load_error(replace(kBase, "[sim]", "[simx]")).find("sim"), std::string::npos);
  // Wrong type and unknown values are errors, not silent defaults.
  EXPECT_NE(load_error(replace(kBase, "mass_kg = 2.0", "mass_kg = \"two\"")).find("body.mass_kg"),
            std::string::npos);
  EXPECT_NE(load_error(replace(kBase, "\"rk4\"", "\"rk5\"")).find("sim.integrator"), std::string::npos);
  EXPECT_NE(load_error(std::string(kBase) + "[[forces]]\ntype = \"warp_drive\"\n").find("forces[0].type"),
            std::string::npos);
}

TEST(Scenario, LoadsExampleScenario) {
  const auto dir = std::filesystem::path(FLIGHTSIM_SOURCE_DIR) / "scenarios";
  const auto s = load_rigid_body_scenario(dir / "rb_linear_drag.toml");
  EXPECT_EQ(s.run.integ, Integrator::RK4);
  EXPECT_EQ(s.mp.mass, 2.0);
  ASSERT_EQ(s.forces.size(), 1u);
  EXPECT_EQ(s.forces[0](0.0, s.x0).F_N, -0.5 * s.x0.segment<3>(3));  // c = 0.5 N s/m
  EXPECT_EQ(s.x0.segment<3>(3), Vec3(10, -5, -20));
  EXPECT_FALSE(s.source_hash.empty());
  // The other shipped examples load too.
  EXPECT_EQ(load_rigid_body_scenario(dir / "rb_torque_free_asym.toml").x0.segment<3>(10),
            Vec3(0.5, 0.0, 1.0));
  EXPECT_EQ(load_rigid_body_scenario(dir / "rb_dzhanibekov.toml").x0.segment<3>(10), Vec3(1e-3, 1.0, 1e-3));
}
