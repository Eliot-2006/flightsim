"""Shared helpers for the JSBSim spike (throwaway)."""
import os, jsbsim
ROOT = os.path.dirname(os.path.abspath(__file__))
PLANET = os.path.join(ROOT, "planet", "zero_rot.xml")

def make(planet=False, gravity_model=None, dt=1/120, root=ROOT, model="toy"):
    fdm = jsbsim.FGFDMExec(root)
    fdm.set_debug_level(0)
    if planet:
        assert fdm.load_planet(PLANET, False)
    assert fdm.load_model(model)
    fdm.set_dt(dt)
    if gravity_model is not None:
        fdm["simulation/gravity-model"] = gravity_model
    return fdm

def init(fdm, vt=500.0, h=5000.0, lat=0.0, lon=0.0, psi=0.0, gamma=0.0):
    fdm["ic/vt-fps"] = vt; fdm["ic/h-sl-ft"] = h
    fdm["ic/lat-gc-deg"] = lat; fdm["ic/long-gc-deg"] = lon
    fdm["ic/psi-true-deg"] = psi; fdm["ic/gamma-deg"] = gamma
    fdm["ic/alpha-deg"] = 2.0; fdm["ic/beta-deg"] = 0.0
    fdm["propulsion/set-running"] = -1
    assert fdm.run_ic()
