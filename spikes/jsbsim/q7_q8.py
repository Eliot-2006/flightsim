"""Q7 stock f16 surface drive, Q8 inertia sign with/without negated_crossproduct_inertia."""
import jsbsim, os, re, shutil
import common
from common import *
print("jsbsim", jsbsim.__version__)
# ---- Q7 (stock aircraft dir shipped in the wheel)
fdm = jsbsim.FGFDMExec(None); fdm.set_debug_level(0)
print("Q7 default root:", fdm.get_root_dir())
fdm.load_model("f16"); fdm.set_dt(1/120)
for k, v in {"ic/vt-fps": 500, "ic/h-sl-ft": 5000, "ic/alpha-deg": 2, "ic/beta-deg": 0, "ic/gamma-deg":0}.items(): fdm[k] = v
fdm["propulsion/set-running"] = -1; fdm.run_ic()
print(" stock f16 weight lbs", fdm["inertia/weight-lbs"], "ixz", fdm["inertia/ixz-slugs_ft2"])
def probe(prop, val, steps=1):
    f = jsbsim.FGFDMExec(None); f.set_debug_level(0); f.load_model("f16"); f.set_dt(1/120)
    for k, v in {"ic/vt-fps": 500, "ic/h-sl-ft": 5000, "ic/alpha-deg": 2, "ic/beta-deg": 0, "ic/gamma-deg":0}.items(): f[k] = v
    f["propulsion/set-running"] = -1; f.run_ic()
    f[prop] = val
    for _ in range(steps): f.run()
    return f["fcs/elevator-pos-rad"], f["fcs/elevator-pos-norm"]
for prop, val in (("fcs/elevator-pos-rad", 0.2), ("fcs/elevator-pos-norm", 0.5), ("fcs/elevator-cmd-norm", 0.5)):
    r1 = probe(prop, val, 1); r120 = probe(prop, val, 120)
    print(f" write {prop}={val}: after 1 step pos_rad,pos_norm={r1[0]:+.4f},{r1[1]:+.4f}; after 1 s={r120[0]:+.4f},{r120[1]:+.4f}")
# is the property writable at all / does FCS recompute? (in ref f16.xml, elevator-pos-rad is the output of aerosurface_scale 'fcs/elevator-position')
# ---- Q8: same toy, default negated_crossproduct_inertia
src = open(os.path.join(ROOT, "aircraft/toy/toy.xml")).read()
os.makedirs(os.path.join(ROOT, "aircraft/toy_neg"), exist_ok=True)
open(os.path.join(ROOT, "aircraft/toy_neg/toy_neg.xml"), "w").write(
    src.replace(' negated_crossproduct_inertia="false"', "").replace('name="toy"', 'name="toy_neg"'))
for m in ("toy", "toy_neg"):
    f = make(model=m); init(f)
    print(f"Q8 {m}: file ixz=+982 -> inertia/ixz-slugs_ft2 = {f['inertia/ixz-slugs_ft2']}")
# physical check: roll-rate IC, free response of r-dot with all aero moments removed would be ideal; instead compare initial rdot
for m in ("toy", "toy_neg"):
    f = make(model=m); init(f)
    f["ic/p-rad_sec"] = 0.5; f.run_ic(); f.run()
    print(f"   {m}: initial rdot with p=0.5 rad/s (tensor Ixz sign differs): rdot={f['accelerations/rdot-rad_sec2']:+.6f} qdot={f['accelerations/qdot-rad_sec2']:+.6f}")
