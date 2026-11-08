"""Q2 load_planet, Q3 gravity-model, Q4 integrator properties."""
import jsbsim, math, common
from common import *
print("jsbsim", jsbsim.__version__)
print("Q2 API: FGFDMExec.load_planet(PlanetPath: str, useAircraftPath: bool) -> bool ; file root <planet>")
for planet in (False, True):
    fdm = make(planet=planet); init(fdm, vt=500, h=5000, lat=45.0)
    print(f" planet={planet}: gravity-model default={fdm['simulation/gravity-model']}  r_veh={fdm['position/radius-to-vehicle-ft']:.3f} ft"
          f"  |v_eci|={fdm['velocities/eci-velocity-mag-fps']:.4f}  vt={fdm['velocities/vt-fps']:.4f}"
          f"  omega_planet(eci-vs-ecef proxy) diff={fdm['velocities/eci-velocity-mag-fps']-fdm['velocities/vg-fps']:.4f} fps")
# bad-file behaviour: loading planet AFTER model
fdm = jsbsim.FGFDMExec(ROOT); fdm.load_model("toy")
try: print(" load_planet after load_model ->", fdm.load_planet(PLANET, False))
except Exception as e: print(" load_planet after load_model raised:", e)
# Q3
fdm = make(planet=True); init(fdm)
for gm in (0, 1):
    fdm["simulation/gravity-model"] = gm; fdm.run()
    print(f"Q3 set simulation/gravity-model={gm} -> read {fdm['simulation/gravity-model']}; g at lat0,h=5000ft: {fdm['accelerations/gravity-ft_sec2'] if 'accelerations/gravity-ft_sec2' in fdm.query_property_catalog('') .split('\n') else 'n/a'}")
print(" gravity props:", [l for l in fdm.query_property_catalog("grav").splitlines()])
# Q4
names = ["simulation/integrator/rate/rotational","simulation/integrator/rate/translational",
         "simulation/integrator/position/rotational","simulation/integrator/position/translational"]
fdm = make(planet=True); init(fdm)
print("Q4 defaults:", {n.split('integrator/')[1]: fdm[n] for n in names})
for n in names:
    ok = {}
    for v in range(-1, 11):
        f = make(planet=True); init(f)
        try: f[n] = v; ok[v] = f[n]
        except Exception as e: ok[v] = "ERR"
    print(f"   {n}: set->read {ok}")
# stability / accuracy of every value on the toy trim, 20 s; ref = AB3/AB2 @ 1/2400 s? use RK-free: compare each to dt=1/4800 default integrators
def run(vals, dt, T=20.0):
    f = make(planet=True, dt=dt); init(f, vt=500, h=5000)
    f["fcs/throttle-cmd-norm[0]"]=0.18025; f["fcs/pitch-trim-cmd-norm"]=-0.04589
    for n, v in zip(names, vals): f[n] = v
    for _ in range(round(T/dt)): f.run()
    return f["position/h-sl-ft"], f["velocities/vt-fps"], f["attitude/theta-deg"]
ref = run((1,3,1,4), 1/4800)
print("   ref (defaults, dt=1/4800): h,vt,theta =", ref)
for v in (1,2,3,4,5,6,7,8,9):
    for kind, idx in (("translational rate",1),("rotational rate",0)):
        vals=[1,3,1,4]; vals[idx]=v
        try:
            r=run(vals,1/120); d=[a-b for a,b in zip(r,ref)]
            print(f"   {kind}={v}, dt=1/120: dh={d[0]:+.4f} dV={d[1]:+.5f} dtheta={d[2]:+.5f}")
        except Exception as e: print(f"   {kind}={v}: ERR {e}")
