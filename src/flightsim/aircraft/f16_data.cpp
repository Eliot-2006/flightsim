#include "flightsim/aircraft/f16_data.hpp"

#include <set>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>

#include "flightsim/math/units.hpp"

namespace flightsim {

namespace {

using units::deg;
using units::ft;
using units::lbf;
using units::slug;

// Reads leaf tables ("section.name", each carrying a `source` key) and remembers which were consumed so that
// finish() can reject tables the loader does not know about.
class Reader {
 public:
  Reader(const std::filesystem::path& p) : file_(p.string()) {
    try {
      root_ = toml::parse_file(file_);
    } catch (const toml::parse_error& e) {
      throw std::runtime_error(file_ + ": " + std::string(e.description()));
    }
    if (root_["schema_version"].value<int>() != 1) fail("schema_version", "expected 1");
  }

  [[noreturn]] void fail(const std::string& path, const std::string& msg) const {
    throw std::runtime_error(file_ + ": " + path + ": " + msg);
  }

  const toml::table& leaf(const std::string& path) {
    const toml::table* t = root_.at_path(path).as_table();
    if (!t) fail(path, "missing table");
    if (!t->contains("source")) fail(path, "table has no `source` citation");
    used_.insert(path);
    return *t;
  }

  double num(const toml::table& t, const std::string& path, const std::string& key) const {
    const auto v = t[key].value<double>();
    if (!v) fail(path + "." + key, "missing or not a number");
    return *v;
  }
  std::vector<double> nums(const toml::table& t, const std::string& path, const std::string& key) const {
    const toml::array* a = t[key].as_array();
    if (!a) fail(path + "." + key, "missing or not an array");
    std::vector<double> out;
    for (const auto& e : *a) {
      const auto v = e.value<double>();
      if (!v) fail(path + "." + key, "array element is not a number");
      out.push_back(*v);
    }
    return out;
  }
  std::vector<std::string> strs(const toml::table& t, const std::string& path, const std::string& key) const {
    const toml::array* a = t[key].as_array();
    if (!a) fail(path + "." + key, "missing or not an array");
    std::vector<std::string> out;
    for (const auto& e : *a) {
      const auto v = e.value<std::string>();
      if (!v) fail(path + "." + key, "array element is not a string");
      out.push_back(*v);
    }
    return out;
  }
  void expect_str(const toml::table& t, const std::string& path, const std::string& key,
                  const std::string& want) const {
    const auto v = t[key].value<std::string>();
    if (!v || *v != want)
      fail(path + "." + key, "expected \"" + want + "\", found " + (v ? "\"" + *v + "\"" : "nothing"));
  }
  void expect_strs(const toml::table& t, const std::string& path, const std::string& key,
                   const std::vector<std::string>& want) const {
    if (strs(t, path, key) != want) fail(path + "." + key, "unexpected axis/unit labels");
  }

  // value `key` with its unit label in `<key>_units` (or `units_key`), which must equal `want`; returns
  // value*factor.
  double scalar(const std::string& path, const std::string& key, const std::string& want, double factor,
                const std::string& units_key = "") {
    const toml::table& t = leaf(path);
    expect_str(t, path, units_key.empty() ? key + "_units" : units_key, want);
    return num(t, path, key) * factor;
  }
  // [path] value + units
  double value(const std::string& path, const std::string& want, double factor) {
    const toml::table& t = leaf(path);
    expect_str(t, path, "units", want);
    return num(t, path, "value") * factor;
  }
  // sibling key of an already-open leaf table
  double key(const std::string& path, const std::string& k, const std::string& want, double factor,
             const std::string& units_key = "") {
    const toml::table& t = leaf(path);
    expect_str(t, path, units_key.empty() ? k + "_units" : units_key, want);
    return num(t, path, k) * factor;
  }
  std::vector<double> pair(const std::string& path, const std::string& k, const std::string& units_key,
                           const std::string& want, double factor) {
    const toml::table& t = leaf(path);
    expect_str(t, path, units_key, want);
    auto v = nums(t, path, k);
    for (double& x : v) x *= factor;
    return v;
  }

  // Rows follow axes[0], columns axes[1]. Breakpoints are scaled by f0/f1, values by fv.
  Table2D table2(const std::string& path, const std::string& ax0, const std::string& ax1,
                 const std::string& u0, const std::string& u1, double f0, double f1,
                 const std::string& vunits = "1", double fv = 1.0) {
    const toml::table& t = leaf(path);
    expect_strs(t, path, "axes", {ax0, ax1});
    expect_strs(t, path, "axis_units", {u0, u1});
    expect_str(t, path, "units", vunits);
    auto x = nums(t, path, ax0);
    auto y = nums(t, path, ax1);
    for (double& v : x) v *= f0;
    for (double& v : y) v *= f1;
    const Eigen::MatrixXd m = matrix(t, path, x.size(), y.size(), fv);
    try {
      return Table2D(path, std::move(x), std::move(y), m);
    } catch (const std::invalid_argument& e) {
      fail(path, e.what());
    }
  }

  Table1D table1(const std::string& path) {
    const toml::table& t = leaf(path);
    expect_strs(t, path, "axes", {"alpha"});
    expect_strs(t, path, "axis_units", {"deg"});
    expect_str(t, path, "units", "1");
    auto x = nums(t, path, "alpha");
    for (double& v : x) v *= deg;
    try {
      return Table1D(path, std::move(x), nums(t, path, "values"));
    } catch (const std::invalid_argument& e) {
      fail(path, e.what());
    }
  }

  // [aero.damping]: columns, in the stated order, become one Table1D(alpha) each.
  std::array<Table1D, 9> damping() {
    const std::string path = "aero.damping";
    const toml::table& t = leaf(path);
    const std::vector<std::string> names = {"CXq", "CYr", "CYp", "CZq", "Clr", "Clp", "Cmq", "Cnr", "Cnp"};
    expect_strs(t, path, "axes", {"alpha", "derivative"});
    expect_strs(t, path, "axis_units", {"deg", "1"});
    expect_strs(t, path, "derivative", names);
    expect_str(t, path, "units", "1");
    auto x = nums(t, path, "alpha");
    for (double& v : x) v *= deg;
    const Eigen::MatrixXd m = matrix(t, path, x.size(), names.size(), 1.0);
    auto col = [&](int j) {
      const std::vector<double> y(m.col(j).data(), m.col(j).data() + m.rows());
      try {
        return Table1D(path + "." + names[static_cast<std::size_t>(j)], x, y);
      } catch (const std::invalid_argument& e) {
        fail(path, e.what());
      }
    };
    return {col(0), col(1), col(2), col(3), col(4), col(5), col(6), col(7), col(8)};
  }

  F16Actuator actuator(const std::string& path) {
    const auto lim = pair(path, "position_limits", "position_units", "deg", deg);
    if (lim.size() != 2) fail(path + ".position_limits", "expected 2 values");
    return {lim[0], lim[1], key(path, "rate_limit", "deg/s", deg, "rate_units"),
            key(path, "lag_rate", "1/s", 1.0)};
  }

  // Rejects any leaf table in the file that no loader call consumed.
  void finish() const { walk(root_, ""); }

 private:
  Eigen::MatrixXd matrix(const toml::table& t, const std::string& path, std::size_t rows, std::size_t cols,
                         double fv) const {
    const toml::array* a = t["values"].as_array();
    if (!a) fail(path + ".values", "missing or not an array");
    if (a->size() != rows) fail(path + ".values", "row count does not match first axis");
    Eigen::MatrixXd m(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
    for (std::size_t i = 0; i < rows; ++i) {
      const toml::array* r = a->get(i)->as_array();
      if (!r || r->size() != cols)
        fail(path + ".values", "row " + std::to_string(i) + " length does not match second axis");
      for (std::size_t j = 0; j < cols; ++j) {
        const auto v = r->get(j)->value<double>();
        if (!v) fail(path + ".values", "non-numeric entry");
        m(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = *v * fv;
      }
    }
    return m;
  }

  void walk(const toml::table& t, const std::string& prefix) const {
    for (const auto& [k, node] : t) {
      const toml::table* sub = node.as_table();
      if (!sub) continue;
      const std::string path = prefix.empty() ? std::string(k.str()) : prefix + "." + std::string(k.str());
      if (sub->contains("source")) {
        if (!used_.count(path)) fail(path, "unknown table (not part of the F16Data schema)");
      } else {
        walk(*sub, path);
      }
    }
  }

  std::string file_;
  toml::table root_;
  std::set<std::string> used_;
};

std::array<double, 3> levels(Reader& r) {
  const auto v = r.pair("engine.thrust_power_levels", "values", "units", "percent", 1.0);
  if (v.size() != 3) r.fail("engine.thrust_power_levels.values", "expected 3 values");
  return {v[0], v[1], v[2]};
}

}  // namespace

F16Data load_f16_data(const std::filesystem::path& toml_path) {
  Reader r(toml_path);
  const double i2 = slug * ft * ft;  // slug*ft^2 -> kg*m^2
  const std::string I = "slug*ft^2";
  auto dmp = r.damping();  // columns -> 9 Table1D, in the TOML's stated derivative order

  F16Data d{
      .weight_nominal = r.value("mass_geometry.W", "lbf", lbf),
      // S&L: mass from the reciprocal-mass constant, not W/g.
      .mass = slug / r.value("mass_geometry.reciprocal_mass", "1/slug", 1.0),
      .Ixx = r.value("mass_geometry.Ixx", I, i2),
      .Iyy = r.value("mass_geometry.Iyy", I, i2),
      .Izz = r.value("mass_geometry.Izz", I, i2),
      .Ixz = r.value("mass_geometry.Ixz", I, i2),
      .S = r.value("mass_geometry.S", "ft^2", ft * ft),
      .b = r.value("mass_geometry.b", "ft", ft),
      .cbar = r.value("mass_geometry.cbar", "ft", ft),
      .xcg_ref = r.value("mass_geometry.xcg_ref", "cbar", 1.0),
      .hx = r.value("mass_geometry.hx", "slug*ft^2/s", i2),
      .g_sl = r.value("environment.g", "ft/s^2", ft),
      .rtod = r.value("environment.rtod", "deg/rad", 1.0),
      .adc = {.rho0 = r.value("environment.rho0", "slug/ft^3", 1.0),
              .temperature0 = r.value("environment.temperature0", "degR", 1.0),
              .temperature_high = r.value("environment.temperature_high", "degR", 1.0),
              .break_altitude = r.value("environment.temperature_break_altitude", "ft", 1.0),
              .altitude_factor = r.value("environment.temperature_altitude_factor", "1/ft", 1.0),
              .density_exponent = r.value("environment.density_exponent", "1", 1.0),
              .gamma = r.value("environment.gamma_air", "1", 1.0),
              .gas_constant = r.value("environment.gas_constant", "ft^2/(s^2*degR)", 1.0)},

      .cx = r.table2("aero.cx", "alpha", "elevator", "deg", "deg", deg, deg),
      .cm = r.table2("aero.cm", "alpha", "elevator", "deg", "deg", deg, deg),
      .cl = r.table2("aero.cl", "alpha", "abs_beta", "deg", "deg", deg, deg),
      .cn = r.table2("aero.cn", "alpha", "abs_beta", "deg", "deg", deg, deg),
      .dlda = r.table2("aero.dlda", "alpha", "beta", "deg", "deg", deg, deg),
      .dldr = r.table2("aero.dldr", "alpha", "beta", "deg", "deg", deg, deg),
      .dnda = r.table2("aero.dnda", "alpha", "beta", "deg", "deg", deg, deg),
      .dndr = r.table2("aero.dndr", "alpha", "beta", "deg", "deg", deg, deg),
      .cz = r.table1("aero.cz"),
      .cz_beta_scale = r.scalar("aero.cz", "beta_scale", "deg", deg),
      .cz_elevator_gain = r.scalar("aero.cz", "elevator_gain", "1", 1.0),
      .cz_elevator_scale = r.scalar("aero.cz", "elevator_scale", "deg", deg),
      .cy_beta_gain = r.scalar("aero.cy", "beta_gain", "1/deg", 1.0 / deg),
      .cy_aileron_gain = r.scalar("aero.cy", "aileron_gain", "1", 1.0),
      .cy_rudder_gain = r.scalar("aero.cy", "rudder_gain", "1", 1.0),
      .aileron_norm = r.scalar("aero.control_normalization", "aileron", "deg", deg, "units"),
      .rudder_norm = r.scalar("aero.control_normalization", "rudder", "deg", deg, "units"),

      .cxq = std::move(dmp[0]),
      .cyr = std::move(dmp[1]),
      .cyp = std::move(dmp[2]),
      .czq = std::move(dmp[3]),
      .clr = std::move(dmp[4]),
      .clp = std::move(dmp[5]),
      .cmq = std::move(dmp[6]),
      .cnr = std::move(dmp[7]),
      .cnp = std::move(dmp[8]),

      // thrust: altitude ft -> m, Mach unchanged, lbf -> N
      .thrust_idle = r.table2("engine.idle", "altitude", "mach", "ft", "1", ft, 1.0, "lbf", lbf),
      .thrust_mil = r.table2("engine.military", "altitude", "mach", "ft", "1", ft, 1.0, "lbf", lbf),
      .thrust_max = r.table2("engine.maximum", "altitude", "mach", "ft", "1", ft, 1.0, "lbf", lbf),
      .power_levels = levels(r),

      .gear_breakpoint = r.scalar("engine.throttle_gearing", "breakpoint", "1", 1.0),
      .gear_low_slope = r.scalar("engine.throttle_gearing", "low_slope", "percent", 1.0, "slope_units"),
      .gear_high_slope = r.scalar("engine.throttle_gearing", "high_slope", "percent", 1.0, "slope_units"),
      .gear_high_intercept =
          r.scalar("engine.throttle_gearing", "high_intercept", "percent", 1.0, "intercept_units"),
      .rtau_dp_low = r.pair("engine.rtau", "delta_power", "delta_power_units", "percent", 1.0).at(0),
      .rtau_dp_high = r.pair("engine.rtau", "delta_power", "delta_power_units", "percent", 1.0).at(1),
      .rtau_low = r.pair("engine.rtau", "values", "units", "1/s", 1.0).at(0),
      .rtau_high = r.pair("engine.rtau", "values", "units", "1/s", 1.0).at(1),
      .rtau_mid_intercept = r.scalar("engine.rtau", "middle_intercept", "1/s", 1.0),
      .rtau_mid_slope = r.scalar("engine.rtau", "middle_slope", "1/(s*percent)", 1.0),
      .power_threshold = r.scalar("engine.power_lag", "power_threshold", "percent", 1.0, "power_units"),
      .power_up_target = r.scalar("engine.power_lag", "transition_up_target", "percent", 1.0, "power_units"),
      .power_down_target =
          r.scalar("engine.power_lag", "transition_down_target", "percent", 1.0, "power_units"),
      .power_high_rate = r.scalar("engine.power_lag", "high_power_rate", "1/s", 1.0),

      .elevator = r.actuator("actuators.elevator"),
      .aileron = r.actuator("actuators.aileron"),
      .rudder = r.actuator("actuators.rudder"),
      .throttle_min = r.pair("actuators.throttle", "position_limits", "units", "1", 1.0).at(0),
      .throttle_max = r.pair("actuators.throttle", "position_limits", "units", "1", 1.0).at(1),
  };
  r.leaf("actuators.aerobench");  // documented as validated-present only
  r.finish();
  return d;
}

}  // namespace flightsim
