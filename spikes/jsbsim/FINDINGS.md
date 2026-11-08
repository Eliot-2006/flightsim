# JSBSim feasibility spike (Task 4.0) — findings

Throwaway. Env: `uv` project in this dir, Python 3.12.13 (uv-managed; system `/usr/local/bin/python3*` are x86 and unusable on this arm64 Mac), **jsbsim 1.3.3** (PyPI wheel, pinned by `uv.lock`), numpy 2.5.4.
Run: `cd spikes/jsbsim && .venv/bin/python <script>`; captured outputs are `out_*.txt`. Files: `aircraft/toy/toy.xml` (toy aircraft), `planet/zero_rot.xml`, `engine/*.xml`, `common.py`.

## Answers

| # | Question | Yes/No | Exact API / property | Evidence |
|---|---|---|---|---|
| 1 | `FGFDMExec(root_dir)` loads custom aircraft from custom dir | **Yes** | `jsbsim.FGFDMExec(root)`; `fdm.load_model("toy")` reads `<root>/aircraft/toy/toy.xml`; engines from `<root>/engine/` | `q1_5_6_8.py` → `out_q1_5_6_8.txt` |
| 2 | Custom planet with `rotation_rate=0` loads from Python | **Yes** | `fdm.load_planet(path: str, useAircraftPath: bool) -> bool` (v1.3.3; C++ `FGFDMExec::LoadPlanet`). File root `<planet name=..>` with `semimajor_axis`, `semiminor_axis`, `rotation_rate unit="RAD/SEC"`, `GM unit="FT3/SEC2"`, `J2`. Verified effect: \|v_eci\|=vt (500.0000) vs 1187.80 with default Earth; radius at lat 45°, 5000 ft = a+5000 | `planet/zero_rot.xml`, `q2_3_4.py` → `out_q2_3_4.txt` |
| 3 | `simulation/gravity-model` settable to standard (0) | **Yes** | `fdm["simulation/gravity-model"] = 0` (default is **1** = WGS84/J2). Std model: g = GM/r² (32.1313 ft/s² at a+5000 ft) | `out_q2_3_4.txt`, `out_earth_effect.txt` |
| 4 | Integrator properties settable | **Yes** (with caveats) | `simulation/integrator/{rate,position}/{rotational,translational}`. Enum: 0 None(freeze), 1 RectEuler, 2 Trapezoidal, 3 AdamsBashforth2, 4 AB3, 5 AB4, 6 Buss1, 7 Buss2, 8 LocalLinearization, 9 AB5. **Defaults: rate/rot=1 (Euler), rate/trans=3 (AB2), position/rot=1 (Euler), position/trans=4 (AB3)**. Set accepts any number; invalid combos throw at `run()`: 6/7/8 are only valid for `position/rotational` (error "Can only use Buss… for rotational position"). No RK4 exists | `out_q2_3_4.txt`; source `FGPropagate.cpp` (v1.3.3) |
| 5 | Minimal custom aircraft with body-axis aero + external force flies and trims | **Yes** | `<axis name="X|Y|Z|ROLL|PITCH|YAW">` functions; `<external_reactions><force frame="BODY">` with `<function>` magnitude; trim: `fdm["simulation/do_simple_trim"]=1` (full trim; raises `TrimFailureError` on failure). Trim result: α=2.747°, throttle 0.18025, pitch trim −0.04589, udot/wdot/qdot ≈ 0/−1.9e-4/4e-16; 10 s after trim Δh=0.57 ft | `out_q1_5_6_8.txt` |
| 6 | Surfaces drivable directly (no FCS) | **Yes**, by design | Toy aero uses a user property (`aero/de-direct-rad`, declared in `<flight_control>`) that no component writes: `fdm["aero/de-direct-rad"]=0.1` → `aero/de-total`=0.1000, q=−0.337 rad/s after 0.5 s. Writing an FCS-*output* property (`fcs/elevator-pos-rad`) is overwritten every frame. Pass-through path: `fcs/elevator-cmd-norm`=0.25 → pos 0.1 rad (range ±0.4) | `out_q1_5_6_8.txt` |
| 7 | Stock `f16`: surfaces settable bypassing FCS? Provenance | **No** | `fcs/elevator-pos-rad` write ignored (+0.2 → +0.0207 after 1 step); `fcs/elevator-pos-norm` write passes one step then FCS scheduler/PID/kinematic (0.3 s traverse) overwrites (→ −0.012 after 1 s). Aero reads `fcs/*-pos-rad` which FCS recomputes. Only route is an edited *copy* of f16.xml. Provenance: header "General Dynamics F-16A", author Erik Hofman, **GPL**, refs: Dash-1 (1987), NASA TP-1538 wind-tunnel data, Richard Murray's Caltech F-16 model, ground-effect papers; "no relation to manufacturer"; not a straight S&L transcription | `out_q7_q8.txt`; `aircraft/f16/f16.xml` in wheel |
| 8 | `negated_crossproduct_inertia="false"` accepted, reflected in `inertia/ixz-slugs_ft2` | **Yes** | File `<ixz>982</ixz>`: `"false"` → property **−982.0**; attribute absent (default) → **+982.0**. Property is the inertia-*tensor* entry; with `"false"` it equals −(file value), i.e. file value is the conventional product of inertia ∫xz dm (positive-Ixz convention, matches S&L). Initial ṙ for p=0.5 rad/s flips sign between the two (−0.0212 vs +0.0230) | `out_q7_q8.txt` |

## Earth-model effect (toy aircraft, 20 s, fixed trim controls, 500 ft/s at 5000 ft, dt=1/240)
Zero-rotation spherical planet + `gravity-model=0` minus JSBSim default (WGS84 rotating, model 1). `earth_effect.py` → `out_earth_effect.txt`.

| lat, heading | Δh (ft) | ΔN / ΔE (ft) | Δθ (deg) | ΔV (ft/s) |
|---|---|---|---|---|
| 0°, N | −7.83 | +3.34 / −0.00 | −0.068 | +0.445 |
| 0°, E | −17.47 | 0 / +7.43 | −0.152 | +0.992 |
| 45°, N | +3.52 | −18.22 / −6.76 | +0.031 | −0.199 |
| 45°, E | −3.29 | +6.82 / −15.46 | −0.029 | +0.186 |

Gravity at end: default 32.183 (0°) – 32.213 (45°) ft/s²; std model 32.131. Differences are 3–17 ft in altitude and up to ~18 ft in ground position over 20 s — non-negligible vs the PRD's 2 % criterion, so the planet+gravity matching is required, as planned. (Δψ −360 in the first row is a wrap, not an effect.) Rotation affects results mainly via Coriolis/ECI start velocity (default |v_eci| = 1187.8 vs 500 ft/s).

## Recommendation
**PRD §8.6 primary approach is feasible; the fallback is not needed.** Custom aircraft dir, zero-rotation spherical planet, standard gravity, body-axis X/Y/Z + moment functions, external body-x thrust force, direct surface drive, JSBSim's own trim, and `negated_crossproduct_inertia="false"` all work in jsbsim 1.3.3 from Python. Two spec amendments: (a) JSBSim has **no RK4** — §8.6 "two rates" convergence must use JSBSim's AB/Euler integrators, so run at a very small dt (e.g. 1/4800) as reference and choose integrator values explicitly; (b) stock-`f16` open-loop (J4) cannot drive surfaces directly — J4 must be trim values only, or use an edited copy of f16.xml (GPL, label as derivative).

## Gotchas
- **Units**: JSBSim is imperial internally (ft, slug, lbf, ft·lbf, deg/rad per property name). Planet file units are honoured via `unit=` attributes. Properties are named with units (`-fps`, `-rad`, `-slugs_ft2`).
- **Structural frame**: x aft-positive, inches; `<location>` for AERORP/CG/force in this frame. `inertia/cg-x-in` echoes the structural value. AERORP≠CG makes JSBSim transfer force moments to the CG (not tested here; toy had AERORP=CG).
- **Integrator defaults** differ per channel (Euler, AB2, Euler, AB3); Euler for rotational rate/position is 1st order, so dt=1/120 gave Δh≈−0.9 ft vs dt=1/4800 reference over 20 s on the toy. Rate values 3,4,5,9 are indistinguishable at that dt; values 6–8 error out for rate and translational position.
- **Table edges clamp**: 1-D table (0→0, 0.1→1, 0.2→4): α=−0.1→0, α=0.3→4.0000, 0.5→4.0000 (`out_q_edge.txt`). Constant extrapolation, unlike S&L's one-cell linear — comparison cases must stay in-domain, as the PRD says.
- **Throttle/trim needs an engine**: `FGTrim` drives throttle via `FGFCS::ThrottleCmd[i]` filled only by engines; without an engine "udot doesn't appear to be trimmable" and trim fails. The toy has a dummy zero-power `<electric_engine>` + direct thruster (`engine/dummy_eng.xml`) and thrust from the external force reading `fcs/throttle-cmd-norm[0]`. Trim also needs `fcs/pitch-trim-cmd-norm` (and roll/yaw trim) wired into the surface chain; an `fcs/elevator-cmd-norm` alone is not what FGTrim varies.
- **Model requirements**: `<ground_reactions/>` element required (load fails without it); user `<property>` declarations must be inside `<flight_control>` (top-level ones are not created); properties that already exist (e.g. `fcs/pitch-trim-cmd-norm`) produce a "already defined" warning.
- `do_simple_trim` value 1 runs a *full* trim (log says "Full Trim"); use other values for other modes — not explored.
- `load_planet` re-initialises models (`InitializeModels`); it returned True even after `load_model` — call it before `load_model` and before IC setup to be safe. `gravity-model=1` with `J2` absent/0 logs "WGS84 gravity model has been set without specifying the J2 gravitational constant."
- Default Earth has rotation, oblateness and latitude-dependent g; initial `ic/vt-fps` is relative to the Earth surface, so the ECI velocity differs by up to ~688 ft/s at lat 45° — another reason to use the zero-rotation planet.
- Python 3.12 via uv-managed install: the machine's `/usr/local` pythons are x86 and break `uv` interpreter discovery; `uv python install 3.12` + `uv sync --managed-python` worked.
