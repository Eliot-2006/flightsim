"""Quantify Earth-model effect on the toy aircraft: 20 s straight flight (fixed trim controls, no feedback)."""
import math, common
from common import *
T, DT = 20.0, 1/240
def run(planet, gm, psi, lat):
    f = make(planet=planet, gravity_model=gm, dt=DT)
    init(f, vt=500, h=5000, lat=lat, psi=psi)
    f["fcs/throttle-cmd-norm[0]"] = 0.18025; f["fcs/pitch-trim-cmd-norm"] = -0.04589
    # warm up with zero-length run so both start from the same IC
    for _ in range(round(T/DT)): f.run()
    return dict(h=f["position/h-sl-ft"], lat=f["position/lat-gc-deg"], lon=f["position/long-gc-deg"],
                th=f["attitude/theta-deg"], psi=f["attitude/psi-deg"], phi=f["attitude/phi-deg"],
                vt=f["velocities/vt-fps"], r=f["position/radius-to-vehicle-ft"], g=f["accelerations/gravity-ft_sec2"])
def diff(a, b, lat):
    rad = a["r"]
    dn = math.radians(a["lat"]-b["lat"])*rad
    de = math.radians(a["lon"]-b["lon"])*rad*math.cos(math.radians(lat))
    return dict(dh_ft=a["h"]-b["h"], dnorth_ft=dn, deast_ft=de, dtheta_deg=a["th"]-b["th"], dpsi_deg=a["psi"]-b["psi"], dphi_deg=a["phi"]-b["phi"], dV_fps=a["vt"]-b["vt"])
for lat in (0.0, 45.0):
  for psi in (0.0, 90.0):
    base = run(False, None, psi, lat)          # JSBSim default: WGS84 rotating, gravity-model default (=1 WGS84 J2)
    zr   = run(True, 0, psi, lat)              # zero-rotation spherical planet + gravity-model 0
    zr1  = run(True, 1, psi, lat)              # zero-rot sphere, gravity-model left at default (1)
    dflt0= run(False, 0, psi, lat)             # default planet, gravity-model 0 only
    print(f"== lat={lat} psi={psi}: T={T}s dt=1/240; default final h={base['h']:.3f} lat={base['lat']:.6f} lon={base['lon']:.6f} g={base['g']:.5f}")
    for name, x in (("zero-rot+gm0 minus default", zr), ("zero-rot + gm1 minus default", zr1), ("gm0 only (default planet) minus default", dflt0)):
        d = diff(x, base, lat)
        print(f"   {name}: " + " ".join(f"{k}={v:+.4f}" for k, v in d.items()) + f"  g_end {x['g']:.5f} vs {base['g']:.5f}")
