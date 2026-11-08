#include "flightsim/sim/output.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <toml++/toml.hpp>

#include "flightsim/math/frames.hpp"
#include "flightsim/version.hpp"

#ifndef FLIGHTSIM_VERSION
#error "FLIGHTSIM_VERSION must be defined by the build"
#endif

namespace flightsim::sim {

namespace {

std::string num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

toml::array arr(const Eigen::Ref<const Eigen::VectorXd>& v) {
  toml::array a;
  for (Eigen::Index i = 0; i < v.size(); ++i) a.push_back(v(i));
  return a;
}

std::ofstream open_out(const std::filesystem::path& p) {
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error(p.string() + ": cannot open for writing");
  return f;
}

}  // namespace

void write_csv(const std::filesystem::path& path, const RigidBodyScenario& s, const RunOutput& out) {
  std::ofstream f = open_out(path);
  f << "t_s,pN_m,pE_m,pD_m,vN_mps,vE_mps,vD_mps,q0,q1,q2,q3,p_radps,q_radps,r_radps,"
       "phi_rad,theta_rad,psi_rad,h_m,V_mps,alpha_rad,beta_rad,E_J,HN_x_kgm2ps,HN_y_kgm2ps,HN_z_kgm2ps,"
       "qnorm\n";
  for (const Sample& smp : out.samples) {
    const Vec<13>& x = smp.x;
    const Vec4 q = rb::quat(x);
    const Vec3 v = rb::vel(x), w = rb::omega(x), p = rb::pos(x);
    const Mat3 C = dcm_from_quat(q);
    const Vec3 euler = euler_from_quat(q);
    const AirData air = air_data_from_body_velocity(C.transpose() * v);
    const double E =
        0.5 * s.mp.mass * v.squaredNorm() + 0.5 * w.dot(s.mp.inertia * w) - s.mp.mass * s.g_N.dot(p);
    const Vec3 H = C * (s.mp.inertia * w + s.mp.h_rotor);
    const double row[] = {smp.t, p.x(),     p.y(),    p.z(), v.x(), v.y(),     v.z(),     q(0),      q(1),
                          q(2),  q(3),      w.x(),    w.y(), w.z(), euler.x(), euler.y(), euler.z(), -p.z(),
                          air.V, air.alpha, air.beta, E,     H.x(), H.y(),     H.z(),     q.norm()};
    bool first = true;
    for (double c : row) {
      f << (first ? "" : ",") << num(c);
      first = false;
    }
    f << '\n';
  }
}

void write_metadata(const std::filesystem::path& path, const RigidBodyScenario& s, const RunOutput& out) {
  toml::table meta;
  meta.insert("flightsim_version", FLIGHTSIM_VERSION);
  meta.insert("git_sha", FLIGHTSIM_GIT_SHA);
  meta.insert("build_type", FLIGHTSIM_BUILD_TYPE);
  meta.insert("compiler", FLIGHTSIM_COMPILER);
  meta.insert("scenario_file", s.source_path);
  meta.insert("scenario_hash", s.source_hash);
  meta.insert("integrator", to_string(s.run.integ));
  meta.insert("dt_s", s.run.dt);
  meta.insert("t_final_s", s.run.t_final);
  meta.insert("output_every", s.run.output_every);
  meta.insert("renormalize", s.run.renormalize);
  meta.insert("stop_at_ground", s.run.stop_at_ground);
  meta.insert("mass_kg", s.mp.mass);
  toml::array inertia;
  for (int i = 0; i < 3; ++i) inertia.push_back(arr(Vec3(s.mp.inertia.row(i))));
  meta.insert("inertia_kgm2", inertia);
  meta.insert("rotor_angular_momentum_kgm2ps", arr(s.mp.h_rotor));
  meta.insert("gravity_NED_mps2", arr(s.g_N));
  meta.insert("x0", arr(s.x0));
  meta.insert("stop_reason", sim::to_string(out.reason));
  meta.insert("t_stop_s", out.t_stop);
  meta.insert("nfev", static_cast<std::int64_t>(out.nfev));
  meta.insert("n_samples", static_cast<std::int64_t>(out.samples.size()));
  if (!out.message.empty()) meta.insert("message", out.message);
  try {  // force parameters live only in the file; embed what was run
    meta.insert("scenario", toml::parse_file(s.source_path));
  } catch (const toml::parse_error&) {
  }
  std::ofstream f = open_out(path);
  f << toml::json_formatter{meta} << '\n';
}

}  // namespace flightsim::sim
