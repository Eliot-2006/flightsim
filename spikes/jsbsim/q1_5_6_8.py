"""Q1 custom dir load, Q5 toy flies+trims, Q6 direct surface drive, Q8 inertia sign."""
import jsbsim, common
from common import *
print("jsbsim", jsbsim.__version__)
# Q1: root_dir = spikes/jsbsim, aircraft in ROOT/aircraft/toy/toy.xml
fdm = make(); init(fdm)
print("Q1 loaded model:", fdm.get_model_name(), "| root:", ROOT)
# Q8: inertia property
print("Q8 inertia/ixz-slugs_ft2 =", fdm["inertia/ixz-slugs_ft2"], "(file value +982, negated_crossproduct_inertia=false)")
print("   ixx,iyy,izz =", fdm["inertia/ixx-slugs_ft2"], fdm["inertia/iyy-slugs_ft2"], fdm["inertia/izz-slugs_ft2"])
print("   weight lbs =", fdm["inertia/weight-lbs"], " cg-x-in (struct) =", fdm["inertia/cg-x-in"])
# Q5: trim (longitudinal)
fdm["fcs/throttle-cmd-norm[0]"] = 0.2
fdm["simulation/do_simple_trim"] = 1  # value 1 = full trim (FGTrim tFull); 2 would be longitudinal etc.
print("Q5 trim: alpha_deg=%.4f theta_deg=%.4f pitch_trim=%.5f elev_pos_rad=%.5f throttle=%.5f"
  % (fdm["aero/alpha-deg"], fdm["attitude/theta-deg"], fdm["fcs/pitch-trim-cmd-norm"], fdm["fcs/elevator-pos-rad"], fdm["fcs/throttle-cmd-norm[0]"]))
print("   udot,wdot,qdot = %.2e %.2e %.2e" % (fdm["accelerations/udot-ft_sec2"], fdm["accelerations/wdot-ft_sec2"], fdm["accelerations/qdot-rad_sec2"]))
h0=fdm["position/h-sl-ft"]
for _ in range(120*10): fdm.run()
print("   10 s after trim: dh=%.3f ft, V=%.3f" % (fdm["position/h-sl-ft"]-h0, fdm["velocities/vt-fps"]))
# Q6: direct drive, bypass FCS
fdm = make(); init(fdm)
fdm["fcs/elevator-cmd-norm"]=0.0
fdm["aero/de-direct-rad"] = 0.1
fdm.run(); 
print("Q6 direct: aero/de-total=%.4f (fcs pos=%.4f, direct=0.1)" % (fdm["aero/de-total"], fdm["fcs/elevator-pos-rad"]))
q0=fdm["velocities/q-rad_sec"]
for _ in range(60): fdm.run()
print("   q after 0.5 s with direct de=+0.1: %.5f rad/s (negative = nose-down, Cm de coeff -1.0)" % fdm["velocities/q-rad_sec"])
# writing the FCS-output property directly
fdm = make(); init(fdm); fdm["fcs/elevator-cmd-norm"]=0.0
fdm["fcs/elevator-pos-rad"]=0.1; fdm.run()
print("   write fcs/elevator-pos-rad=0.1 then run: reads back %.4f (component overwrites output -> %s)" % (fdm["fcs/elevator-pos-rad"], "overwritten" if abs(fdm["fcs/elevator-pos-rad"]-0.1)>1e-9 else "held"))
# cmd path
fdm = make(); init(fdm); fdm["fcs/elevator-cmd-norm"]=0.25; fdm.run()
print("   cmd 0.25 -> pos %.4f rad (range +-0.4)" % fdm["fcs/elevator-pos-rad"])
