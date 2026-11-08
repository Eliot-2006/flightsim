#include "flightsim/sim/scenario.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <toml++/toml.hpp>

#include "flightsim/math/units.hpp"

namespace flightsim::sim {

namespace {

// Typed access to a TOML table that remembers which keys were read, so finish() can reject unknown ones.
class Section {
 public:
  Section(const toml::table& t, std::string file, std::string path)
      : t_(t), file_(std::move(file)), path_(std::move(path)) {}

  [[noreturn]] void fail(const std::string& key, const std::string& msg) const {
    throw std::runtime_error(file_ + ": " + (path_.empty() ? key : path_ + "." + key) + ": " + msg);
  }

  const toml::node* find(const std::string& key) {
    used_.push_back(key);
    return t_.get(key);
  }
  bool has(const std::string& key) const { return t_.contains(key); }

  double number(const std::string& key) {
    const toml::node* n = find(key);
    if (!n) fail(key, "missing required key");
    return as_number(key, n);
  }
  double number(const std::string& key, double dflt) { return has(key) ? number(key) : dflt; }

  std::string string(const std::string& key) {
    const toml::node* n = find(key);
    if (!n) fail(key, "missing required key");
    auto v = n->value<std::string>();
    if (!v) fail(key, "expected a string");
    return *v;
  }
  bool boolean(const std::string& key, bool dflt) {
    if (!has(key)) return dflt;
    auto v = find(key)->value<bool>();
    if (!v) fail(key, "expected true or false");
    return *v;
  }
  long integer(const std::string& key, long dflt) {
    if (!has(key)) return dflt;
    auto v = find(key)->value<std::int64_t>();
    if (!v) fail(key, "expected an integer");
    return static_cast<long>(*v);
  }
  Vec3 vec3(const std::string& key) {
    const toml::node* n = find(key);
    if (!n) fail(key, "missing required key");
    return to_vec3(key, n);
  }
  Vec3 vec3(const std::string& key, const Vec3& dflt) { return has(key) ? vec3(key) : dflt; }
  Mat3 mat3(const std::string& key) {
    const toml::node* n = find(key);
    if (!n) fail(key, "missing required key");
    const toml::array* rows = n->as_array();
    if (!rows || rows->size() != 3) fail(key, "expected a 3x3 array of arrays");
    Mat3 M;
    for (std::size_t i = 0; i < 3; ++i)
      M.row(static_cast<Eigen::Index>(i)) = to_vec3(key, rows->get(i)).transpose();
    return M;
  }
  Section table(const std::string& key) {
    const toml::node* n = find(key);
    if (!n) fail(key, "missing required table");
    const toml::table* sub = n->as_table();
    if (!sub) fail(key, "expected a table");
    return Section(*sub, file_, path_.empty() ? key : path_ + "." + key);
  }
  Section element(const toml::table& t, std::string path) const { return Section(t, file_, std::move(path)); }

  // Rejects any key that was never read.
  void finish() const {
    for (const auto& [key, value] : t_) {
      bool seen = false;
      for (const auto& u : used_) seen = seen || u == key.str();
      if (!seen) fail(std::string(key.str()), "unknown key");
    }
  }

 private:
  double as_number(const std::string& key, const toml::node* n) const {
    auto v = n->value<double>();
    if (!v || !std::isfinite(*v)) fail(key, "expected a finite number");
    return *v;
  }
  Vec3 to_vec3(const std::string& key, const toml::node* n) const {
    const toml::array* a = n ? n->as_array() : nullptr;
    if (!a || a->size() != 3) fail(key, "expected an array of 3 numbers");
    return Vec3(as_number(key, a->get(0)), as_number(key, a->get(1)), as_number(key, a->get(2)));
  }

  const toml::table& t_;
  std::string file_, path_;
  std::vector<std::string> used_;
};

std::string fnv1a64_hex(const std::string& bytes) {
  std::uint64_t h = 1469598103934665603ull;  // ponytail: non-cryptographic identity hash; swap for SHA-256 if
                                             // tamper-evidence is ever needed
  for (char c : bytes) h = (h ^ static_cast<unsigned char>(c)) * 1099511628211ull;
  char buf[32];
  std::snprintf(buf, sizeof buf, "fnv1a64:%016llx", static_cast<unsigned long long>(h));
  return buf;
}

ForceModel parse_force(Section& f) {
  const std::string type = f.string("type");
  ForceModel model;
  if (type == "linear_drag") {
    model = linear_drag(f.number("c_kgps"));
  } else if (type == "quadratic_drag") {
    model = quadratic_drag(f.number("k_kgpm"));
  } else if (type == "spring") {
    model = spring(f.vec3("k_Npm"));
  } else if (type == "body_force") {
    model = body_force(f.vec3("force_N"));
  } else if (type == "body_torque") {
    model = body_torque(f.vec3("torque_Nm"));
  } else {
    f.fail("type", "unknown force type '" + type +
                       "' (linear_drag | quadratic_drag | spring | body_force | body_torque)");
  }
  f.finish();
  return model;
}

}  // namespace

RigidBodyScenario load_rigid_body_scenario(const std::filesystem::path& path) {
  const std::string file = path.string();
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(file + ": cannot open scenario file");
  std::stringstream buf;
  buf << in.rdbuf();
  const std::string text = buf.str();

  toml::table root;
  try {
    root = toml::parse(text, file);
  } catch (const toml::parse_error& e) {
    throw std::runtime_error(file + ": TOML parse error: " + std::string(e.description()) + " (line " +
                             std::to_string(e.source().begin.line) + ", column " +
                             std::to_string(e.source().begin.column) + ")");
  }

  Section top(root, file, "");
  RigidBodyScenario s;
  s.source_path = file;
  s.source_hash = fnv1a64_hex(text);

  Section sim = top.table("sim");
  s.run.integ = Integrator::RK4;
  try {
    s.run.integ = integrator_from_string(sim.string("integrator"));
  } catch (const std::invalid_argument& e) {
    sim.fail("integrator", e.what());
  }
  s.run.dt = sim.number("dt_s");
  s.run.t_final = sim.number("t_final_s");
  const long every = sim.integer("output_every", 1);
  s.run.renormalize = sim.boolean("renormalize", true);
  s.run.stop_at_ground = sim.boolean("stop_at_ground", false);
  sim.finish();
  if (s.run.dt <= 0.0) sim.fail("dt_s", "must be positive");
  if (s.run.t_final <= 0.0) sim.fail("t_final_s", "must be positive");
  if (every < 1) sim.fail("output_every", "must be >= 1");
  s.run.output_every = static_cast<int>(every);
  const double ratio = s.run.t_final / s.run.dt;
  const long steps = std::lround(ratio);
  if (steps < 1 || std::abs(ratio - static_cast<double>(steps)) > 1e-9 * static_cast<double>(steps))
    sim.fail("dt_s", "dt_s must divide t_final_s into a whole number of steps (t_final_s/dt_s = " +
                         std::to_string(ratio) + ")");
  if (steps % every != 0)
    sim.fail("output_every", "output interval must divide the run: " + std::to_string(steps) +
                                 " steps is not a multiple of " + std::to_string(every));

  Section body = top.table("body");
  s.mp.mass = body.number("mass_kg");
  s.mp.inertia = body.mat3("inertia_kgm2");
  s.mp.h_rotor = body.vec3("rotor_angular_momentum_kgm2ps", Vec3::Zero());
  body.finish();
  if (s.mp.mass <= 0.0) body.fail("mass_kg", "must be positive");
  if ((s.mp.inertia - s.mp.inertia.transpose()).cwiseAbs().maxCoeff() >
      1e-12 * s.mp.inertia.cwiseAbs().maxCoeff())
    body.fail("inertia_kgm2", "must be symmetric");
  const Eigen::SelfAdjointEigenSolver<Mat3> eig(0.5 * (s.mp.inertia + s.mp.inertia.transpose()),
                                                Eigen::EigenvaluesOnly);
  if (!(eig.eigenvalues()(0) > 1e-12 * eig.eigenvalues()(2)))
    body.fail("inertia_kgm2", "must be positive definite (singular or non-physical inertia)");

  Section init = top.table("initial");
  const Vec3 p0 = init.vec3("position_m"), v0 = init.vec3("velocity_mps");
  const Vec3 euler = init.vec3("euler_deg") * units::deg, w0 = init.vec3("omega_radps");
  init.finish();
  s.x0 << p0, v0, quat_from_euler(euler.x(), euler.y(), euler.z()), w0;

  double g = units::g0;
  if (top.has("environment")) {
    Section env = top.table("environment");
    g = env.number("gravity_mps2", units::g0);
    env.finish();
  }
  s.g_N = Vec3(0.0, 0.0, g);

  if (const toml::node* fn = top.find("forces")) {
    const toml::array* arr = fn->as_array();
    if (!arr) top.fail("forces", "expected an array of tables ([[forces]])");
    for (std::size_t i = 0; i < arr->size(); ++i) {
      const toml::table* ft = arr->get(i)->as_table();
      if (!ft) top.fail("forces", "expected an array of tables ([[forces]])");
      Section f = top.element(*ft, "forces[" + std::to_string(i) + "]");
      s.forces.push_back(parse_force(f));
    }
  }
  top.finish();
  return s;
}

}  // namespace flightsim::sim
