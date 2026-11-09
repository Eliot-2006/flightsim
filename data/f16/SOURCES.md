# F-16 low-fidelity numeric data provenance

Task 4.1, retrieved 2026-11-08. Only numeric facts were transcribed. No book prose,
figures, downloaded functions, or third-party implementation code are included.
The owner has no copy of Stevens, Lewis & Johnson (S&L); textbook section/table
labels below are the secondary sources' labels, not independently verified book
page citations. This is the reduced S&L model, not the full TP-1538 aerodynamic
model or the Morelli polynomial model.

## Sources and licenses

| ID | Source, pinned version, role | License / independence |
|---|---|---|
| N | Nguyen, L. T.; Ogburn, M. E.; Gilbert, W. P.; Kibler, K. S.; Brown, P. W.; Deal, P. L. (1979), *Simulator Study of Stall/Post-Stall Characteristics of a Fighter Airplane with Relaxed Longitudinal Static Stability*, NASA TP-1538. [NTRS 19800005879](https://ntrs.nasa.gov/citations/19800005879), [PDF](https://ntrs.nasa.gov/api/citations/19800005879/downloads/19800005879.pdf). Primary underlying data. | US federal government report, public domain. Table I, printed p.43 / PDF p.49, was visually checked for mass/geometry/limits; Table VI, printed p.93 / PDF p.99, for thrust. Control signs: symbol definitions, printed pp.5–6 / PDF pp.11–12. |
| A | [AeroBenchVVPython](https://github.com/stanleybak/AeroBenchVVPython/tree/94d44641d56c0881d789208fb8a9deb840a422e3), commit `94d44641d56c0881d789208fb8a9deb840a422e3`, `code/aerobench/lowlevel/*.py`. Primary raw transcription. | [GPL v3 LICENSE](https://github.com/stanleybak/AeroBenchVVPython/blob/94d44641d56c0881d789208fb8a9deb840a422e3/LICENSE). A port of AeroBenchVV; its MATLAB parent is not an independent source. |
| A-M | [AeroBenchVV](https://github.com/pheidlauf/AeroBenchVV/tree/cbb549801955b183eb1b0d5059bd1295e705fa3d), commit `cbb549801955b183eb1b0d5059bd1295e705fa3d`, `src/main/F16_Model/NonlinearModel/*.m`, `src/main/utils/getDefaultSettings.m`, `src/main/F16_Model/{controlledF16,applyCtrlSatLimits}.m`. | [GPL v3 LICENSE.txt](https://github.com/pheidlauf/AeroBenchVV/blob/cbb549801955b183eb1b0d5059bd1295e705fa3d/LICENSE.txt). Same lineage as A. Both MATLAB and Python agree on the disputed raw A aerodynamic entries. |
| B | [Archimedes subsonic F-16 tutorial](https://pinetreelabs.io/archimedes/tutorials/f16/f16_00.html), [source directory](https://github.com/PineTreeLabs/archimedes/tree/5a5b82dc3ee82547c58667f68179136a9d0e2914/docs/source/tutorials/f16), commit `5a5b82dc3ee82547c58667f68179136a9d0e2914`. `aero.py`, `engine.py`, `f16_01.md`, `f16.py`, `config.yaml`, `test_f16.py`, `_static/trim_cases.yaml`. Independent raw B file and published reference values. | [GPL-3.0-or-later](https://github.com/PineTreeLabs/archimedes/blob/5a5b82dc3ee82547c58667f68179136a9d0e2914/pyproject.toml), [LICENSE](https://github.com/PineTreeLabs/archimedes/blob/5a5b82dc3ee82547c58667f68179136a9d0e2914/LICENSE). Part 1 describes a direct Fortran textbook port and lists AeroBench as a similar implementation. It knows about AeroBench: `aero.py` contains disabled alternative Cl numbers labelled as AeroBench modifications. Thus it is a separate reported transcription, but complete independence cannot be established from this declaration alone. C provides a pre-AeroBench check. All B files/pages count as ONE reference lineage. |
| C | University of Minnesota, [DARPA SEC F-16 distribution page](https://dept.aem.umn.edu/~balas/darpa_sec/SEC.Software.html), [F16Simulation.tar.gz](https://dept.aem.umn.edu/~balas/darpa_sec/software/F16Simulation.tar.gz). `F16Sim/lofi_F16_AeroData.c`, `F16_Actuator_Library.mdl`, `nlplant.c`. Additional independent numeric corroboration. | No license declaration was found in the downloaded archive; no code reused. The actuator library records creation in November 2002, predating AeroBench's 2017 copyright. This is an older low/high-fidelity distribution, not an AeroBench descendant. Its **numbers**, not its changed lookup algorithm or whole aircraft model, corroborate disputed reduced tables. |
| J (not B) | [F16Model.jl](https://github.com/isrlab/F16Model.jl/tree/8b719078922b8fa786554c2446ba749c4ca7a447), commit `8b719078922b8fa786554c2446ba749c4ca7a447`. Inspected README, `NonlinearF16Model.jl`, interpolation interface. | [MIT](https://github.com/isrlab/F16Model.jl/blob/8b719078922b8fa786554c2446ba749c4ca7a447/LICENSE). It explicitly credits the UMN package. Current version uses larger NASA HDF5 tables, leading-edge-flap effects and direct thrust input, so it cannot independently check this entire reduced-model/engine schema. Not counted as a complete cross-check. |

Scratch downloads and extraction tools reside outside the repository in
`/tmp/f16src`. Downloaded code was read as text, never imported, built, or run.
Array extraction parsed literal numeric lists only; A and B arrays were read
separately (B was not initialized from A). Transposition/column-major reshape
only aligns axes. Scalar arithmetic labelled “derived” in B translates its
published constants into the same schema; it is not a simulated reference.
NASA W and the audit of AeroBench's lack of actuator dynamics supplement BOTH
files and are not claimed as two independent reproductions.

Download SHA-256 values (the UMN URL is not versioned):

| Artifact | SHA-256 |
|---|---|
| NASA PDF | `aae0ece64474291368c0b4c816d3ab327c6100329e6eb030c2f4545d0913feb3` |
| UMN archive | `c26f0b1c206fcaf4af16e8cafe845227da1ce75900266bd53e7f7b95e1f22349` |
| Rendered B Part 1 | `5b05e6c2f8a43c9d16947254a867eac9c1bd87d3c0846d0da9526146f98400fc` |
| Rendered B Part 3 | `9c841ebb53c41588bf871a49bdcaedbf2566565070ceff34d960df8cc54baa52` |

## Schema and model ownership

`f16_sl.toml` is the resolved source of truth. `crosscheck_sourceB.toml` preserves
raw B values, including disagreements. `schema_version = 1`. Every leaf table
has a `source`; scalar tables use `value` and `units`. Array tables have `axes`,
`axis_units`, named breakpoints and `values`. Rows follow the first axis;
columns the second. Damping's second axis names the nine derivatives in order.
All units are original imperial units; the eventual loader owns SI conversion.
`xcg` and `xcg_ref` are fractions of mean chord, not percentage points.

Inventory: Cx and Cm 12×5; Cl and Cn 12×7; DLDA/DLDR/DNDA/DNDR each 12×7;
Cz 12 entries plus beta/elevator terms; analytic Cy gains; damping 12×9;
idle/military/maximum thrust each 6×6; geometry/inertia, reciprocal mass,
engine rotor momentum, atmosphere constants, throttle gearing, rtau, engine
transition targets, surface limits and supplemental actuator lag/rate limits.
There are no copied implementations. Formula comments describe how to interpret
the numbers, not third-party source code.

S&L table lookup angles are degrees, while state alpha/beta and p/q/r are
radians and radians/s. Body frame is forward/right/down. Positive elevator is
nose-down control (trailing edge down); positive aileron is left-roll control
(right trailing edge down, left up); positive rudder is left-yaw control
(trailing edge left). N's symbol definitions confirm the moment directions.
Cl/Cn store nonnegative beta and are odd in beta. Effectiveness tables use
**signed** beta, not abs(beta). Normalize aileron by 20° and rudder by 30°;
20° is a coefficient normalization, not the 21.5° position limit. Cz's 57.3°
beta divisor and the 57.29578 deg/rad state conversion are distinct rounded
constants. CG transfer uses the *total* Cz/Cy after damping.

Aero domains: alpha [-10,45]°, abs_beta [0,30]° for static Cl/Cn,
signed beta [-30,30]° for effectiveness; elevator [-24,24]° for Cx/Cm.
The no-extrapolation intersection is alpha [-10,45]°, beta [-30,30]°,
elevator [-24,24]°. Elevator's physical ±25° extends 1° beyond its table.
Thrust breakpoints span altitude [0,50000] ft and Mach [0,1]; the reduced
model is intended for subsonic use. The source routines extend the boundary
cell linearly, rather than clamp the output; their arithmetic itself does
not impose a one-cell-distance limit. Caller domain policy belongs to Task 4.2.

## Numeric resolutions

Policy: use TP-1538 for directly identifiable matching values; otherwise use
the reduced-table value supported by the separate B and older C reproductions,
or A and C when B disagrees. A-M is evidence of the A lineage, not an extra vote.
This changes 22 raw A cells in the canonical file. It does not silently turn
B into an AeroBench transcription. NASA's much larger Table III is not a
one-to-one tabulation of the reduced model; the disputed reduced aerodynamic
cells below were resolved by reproduction agreement, not claimed as direct
NASA transcriptions.

All 28 raw A/B differences (including supplemental scalar translations) are
recorded below and in the following subsections. Indices in the machine files
are zero-based. Entries that now agree with B need no comparator exception.

| Cell | A (raw) | B (raw) | Canonical | Evidence |
|---|---:|---:|---:|---|
| Cl(alpha=15°, beta=5°) | -0.022 | -0.019 | -0.019 | B and C AL table |
| Cl(alpha=15°, beta=10°) | -0.041 | -0.034 | -0.034 | B and C AL table |
| Cl(alpha=15°, beta=15°) | -0.054 | -0.044 | -0.044 | B and C AL table |
| Cl(alpha=15°, beta=20°) | -0.06 | -0.046 | -0.046 | B and C AL table |
| Cl(alpha=15°, beta=25°) | -0.063 | -0.046 | -0.046 | B and C AL table |
| Cl(alpha=15°, beta=30°) | -0.068 | -0.047 | -0.047 | B and C AL table |
| Cl(alpha=20°, beta=5°) | -0.022 | -0.02 | -0.02 | B and C AL table |
| Cl(alpha=20°, beta=10°) | -0.045 | -0.04 | -0.04 | B and C AL table |
| Cl(alpha=20°, beta=15°) | -0.057 | -0.05 | -0.05 | B and C AL table |
| Cl(alpha=20°, beta=20°) | -0.069 | -0.059 | -0.059 | B and C AL table |
| Cl(alpha=20°, beta=25°) | -0.081 | -0.068 | -0.068 | B and C AL table |
| Cl(alpha=20°, beta=30°) | -0.089 | -0.074 | -0.074 | B and C AL table |
| Cl(alpha=25°, beta=5°) | -0.021 | -0.02 | -0.02 | B and C AL table |
| Cl(alpha=25°, beta=10°) | -0.04 | -0.037 | -0.037 | B and C AL table |
| Cl(alpha=25°, beta=15°) | -0.054 | -0.049 | -0.049 | B and C AL table |
| Cl(alpha=25°, beta=20°) | -0.067 | -0.061 | -0.061 | B and C AL table |
| Cl(alpha=25°, beta=25°) | -0.079 | -0.071 | -0.071 | B and C AL table |
| Cl(alpha=25°, beta=30°) | -0.088 | -0.079 | -0.079 | B and C AL table |
| Cz(alpha=5°) | -0.415 | -0.416 | -0.416 | B and C A table |
| Cz(alpha=20°) | -1.355 | -1.366 | -1.366 | B and C A table |
| CXq(alpha=-5°) | -0.11 | 0.11 | 0.11 | B and C damping |
| CYr(alpha=40°) | -0.493 | 0.493 | -0.493 | A and C damping |
| Idle thrust(h=10000 ft, M=0.6), lbf | -170 | -710 | -710 | N Table VI: -710; B agrees |

A separate literal-array audit found all 744 canonical aerodynamic cells equal
to C's low-fidelity numeric arrays. C's control-effectiveness routine binds beta
differently; this audit checks numeric arrays, not that algorithm.

C's low-fidelity arrays corroborate B's Cz values -0.416 and -1.366 and all
18 Cl alternatives above. B explicitly distinguishes its textbook Cl rows
from the disabled AeroBench rows. C and B have CXq(-5°)=+0.110; C and A have
CYr(40°)=-0.493. Cmq(-5°)=-0.540 is suspiciously small compared with its
neighbors, but **all three transcriptions agree**, so it is retained; B also
questions this entry. Do not “smooth” it or replace it with -5.40 without new
evidence. These decisions are numerical provenance choices, not claims that
all variants describe identical aircraft behavior.

### Mass and inertia

NASA nominal W is 20500 lbf. A instead multiplies forces by the printed
`rm=0.00157 1/slug`, implying mass 636.9426751592357 slug and effective weight
20490.44585987261 lbf at g=32.17 ft/s². B's Part 1 uses weight 20490.4459 lbf,
so its independently derived reciprocal mass is 0.0015699999969253965 1/slug.
The relative difference is about 1.96e-9, exceeding the required 1e-9 tolerance;
`resolutions.toml` explicitly retains 0.00157. Both files keep NASA's nominal
W separately. **Do not compute check-case mass as W/g**, which would change the
reference implementation. B's `f16.py` comment mentions 25000 lbf as a textbook
weight; this is inconsistent with N and B's own effective weight, so it is
recorded here rather than adopted. C/J use rounded mass 636.94 slug, not an
extra exact match to either reciprocal.

The positive product Ixz=982 slug·ft² maps to **negative** xz/zx entries in the
FRD inertia tensor; B stores that matrix element as -982. Canonical data store
the positive product, not the matrix element. Rotor momentum hx=160 is along
positive body x. A's implementation also has rounded, combined rotational
coefficients (c1…c9); these are not replacements for the published inertias.
The full-inertia EOM derives the coupling from I; the resulting rounding differences
against A1 are investigated in `docs/design.md` (Validation findings), without altering the inertias.

### Engine and atmosphere

The idle-thrust -170 in A is corrected to -710 lbf from visually checked
TP-1538 Table VI (Mach 0.6, altitude 10000 ft); SI in the same table gives
-3158 N, corroborating the magnitude. N Table VI begins at Mach 0.2;
**Mach-zero rows come only from the S&L reproductions A and B**. The other
Mach≥0.2 thrust values were visually checked against Table VI and agree.
Negative idle thrust entries are intentional data, not transcription errors.

Throttle gearing uses low slope 64.94, high slope 217.38 and breakpoint 0.77.
A's high intercept is -117.38 percentage points. B computes -117.3788 from
its slopes/breakpoint to enforce continuity. Canonical retains A's printed
constant; exact values are pinned in `resolutions.toml`. At throttle 0.77 use
the low branch (≤). Power levels 0/50/100 identify idle/military/maximum;
thrust blends linearly on either side of 50. `rtau` is an **inverse** time
constant in 1/s, not a time in seconds: 1 for delta-power ≤25, 0.1 for ≥50,
and 1.9−0.036·delta-power between. For power≥50 use 5 1/s; when crossing the
50 threshold, intermediate targets are 60 on acceleration and 40 on deceleration.
These branch meanings accompany the numerical parameters without copying code.

The simplified atmosphere constants were checked against B's **Part 1 adc**,
which retains density exponent 4.14. B's later `LinearAtmosphere` derives an
exponent from g/R/lapse rate instead, so that class is not an exact S&L adc
cross-check. A changes temperature at altitude ≥35000 ft; B Part 1 uses >35000.
Use A's ≥ boundary in the later S&L atmosphere implementation; the numeric
breakpoint alone does not encode this inequality. A also substitutes altitude
0.01 ft for a negative-altitude thrust query. Do not silently generalize these
implementation quirks into changes to the tabulated data.

### Actuators

A-M `getDefaultSettings.m` and A's `LowLevelController.CtrlLimits` use elevator
±25°, aileron ±21.5°, rudder ±30° and throttle [0,1]. Inspection of MATLAB
`controlledF16.m`, `applyCtrlSatLimits.m`, the Simulink block diagram, and
Python `low_level_controller.py` / `controlled_f16.py` found **instantaneous
position saturation**, with no surface-position states or surface lag/rate
blocks. Thus there are no “lag/rate limits used in AeroBench” to transcribe;
the `actuators.aerobench` flags explicitly record their absence.

For the actuator model, C's `F16_Actuator_Library.mdl` supplies the
S&L-style gain **20.2 1/s** and rates **60/80/120 deg/s** for elevator/aileron/
rudder. B Part 2/config reproduces these rates and positions but uses
**tau=0.0495 s**, equivalent to **20.2020202020202 1/s**. Three exact
`resolutions.toml` entries retain C's 20.2. This is a rounding disagreement,
not an excuse to widen tolerance. NASA Table I independently confirms the
position limits. The data distinguish supplemental S&L/UMN parameters from
what AeroBench actually simulates.

## Published reference values

`reference_values.toml` records all 13 flight-variable derivatives of the S&L
Table 3.5-2 check case, every input/state/control value, plus B's reported
body velocities and (du,dv,dw)=(100.8536,-218.3080,-437.0399) ft/s². The latter
are used to check the NED-velocity state implementation. Values are cited
individually; citation IDs resolve to pinned paths or the rendered tutorial.
The rendered Part 3 page is used for displayed linear-model numbers, rather
than executing its tutorial to generate new numbers.

All published reference outputs found here are from **one independent lineage,
B**, even when two B paths quote them. Every scalar/array accordingly has
`single_source=true`; none is represented as an independently verified book
transcription. A/A-M and C/J supply data, not these printed output fixtures.

Reference qualifications/disagreements:

- **Altitude sign:** B Part 1's expected 12th derivative is hdot=+248.1241
  ft/s. B's Part 3/test uses NED position with pd_dot=-248.1241 ft/s. These
  agree under h=-p_D. B also contains a confusing note describing a DCM typo
  and a sign change (it calls this a yaw-rate issue, though index 11 is altitude).
  Do not negate hdot based on that note. The two stored chart values are
  explicitly labelled; direct FRD→NED kinematics independently confirm them.
- **TR1:** sea-level 502 ft/s level-flight results at xcg=0.35, 0.30 and 0.38
  are transcribed from B's reproduction of Table 3.6-3. In particular at 0.35:
  alpha=0.03691 rad, throttle=0.1385, elevator=-0.7588°. Tiny lateral trim
  values are retained. Table 3.6-2's coordinated turn at xcg=0.35 is a separate
  case. B's YAML has another turn at xcg=0.30; these are different cases, not
  conflicting outputs. No published altitude level-flight trim was located.
- **LN1:** a published matrix at exactly xcg=0.35 was **not located**. B's
  502 ft/s sea-level stability example uses **xcg=0.30**, confirmed by
  `config.yaml` and `test_lon_stability` / `test_lat_stability`. It must not be
  relabelled 0.35. The available longitudinal tutorial A/B matrices (4 states,
  algebraic engine), textbook-quoted lateral matrix (4 states), and textbook-
  quoted eigenvalues are retained under `LN1.level_30`. The exact requested
  0.35 matrix is marked unavailable, not synthesized by finite differences.
- **Eigenvalue disagreement:** textbook values quoted by B2 vs tutorial
  outputs in B3: phugoid -0.0087±0.0739i vs -0.0087±0.0740i; short period
  -1.2039±1.4922i vs -1.2038±1.4920i; Dutch roll -0.4399±3.220i vs
  -0.4362±3.2190i; roll subsidence -3.601 vs -3.6003; spiral both -0.0128.
  Both sets are retained, with their status labelled. No tolerance was chosen
  here to erase these differences. Compare each in the later analysis task.

A1/TR1/LN1 acceptance tolerances belong to their eventual tests (PRD §8.5),
with comments citing printed precision. The transcription comparator's 1e-9
relative tolerance comes directly from this task, **zero absolute tolerance**;
it must not be reused as a tolerance for printed reference outputs.

## Reproduction and checks

Requires Python ≥3.11, standard library plus `tomllib` only. No installation or
model execution is needed:

```sh
python3 python/tools/f16_data_diff.py --self-test
python3 python/tools/f16_data_diff.py data/f16/f16_sl.toml data/f16/crosscheck_sourceB.toml
```

The second command automatically reads `resolutions.toml` beside the canonical
left file; `--resolutions PATH` selects another exact exception list. Each
exception pins path, both values, canonical choice, reason and evidence. Stale,
duplicate, undocumented, changed-value or nonnumeric exceptions fail; array
lengths, missing keys, units and booleans are compared as well. Only `source`
fields are omitted from equality checks because provenance differs. There are
six remaining differences, each printed as `RESOLVED`; any other difference
returns 1. Input/resolution errors return 2. The final line is:

```text
PASS: 6 differences, 6 documented resolutions, 0 unresolved
```

No C++ model, loader, trim solver or linearizer exists on this branch yet, so
these checks establish transcription consistency and tool behavior, not flight-
dynamics correctness. The book and an independent second source for the output
fixtures remain the main limitations for the downstream validation tasks.
