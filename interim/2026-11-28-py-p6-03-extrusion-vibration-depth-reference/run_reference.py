"""Small independent references. Writes a NEW directory; never accepts P6.

Larger --levels/--lengths/--ring-levels/--periods are intended for a user PC.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from elastic_reference import ring, foundation, lame
from vibration_reference import diagnose, exact
from compare_outputs import CONVENTIONS, SCHEMA, FIELDS, validate


def write_json(path, data):
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(data,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write("\n")


def dataset(name, conditions, records, dynamic=False):
    d = {"schema":SCHEMA,"complete":True,"conventions":CONVENTIONS,
         "case":{"name":name,"physical_conditions":conditions,
                 "time_meaning":"physical_seconds" if dynamic else "static_state_label_not_physical_time"},
         "records":records}
    validate(d)
    return d


def record(identifier, xyz, quantity, component, unit, value, t=0., location="physical_point"):
    return dict(zip(FIELDS,[identifier,location,float(t),*map(float,xyz),quantity,component,unit,float(value)]))


def export(out, name, d):
    write_json(out/(name+".json"),d)
    write_json(out/(name+".meta.json"),{k:v for k,v in d.items() if k != "records"})
    with (out/(name+".csv")).open("x",encoding="utf-8",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=FIELDS,lineterminator="\n")
        writer.writeheader()
        writer.writerows(d["records"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    p.add_argument("--levels",nargs="+",type=int,default=[1,2])
    p.add_argument("--lengths",nargs="+",type=float,default=[2.,4.,8.])
    p.add_argument("--widths",nargs="+",type=float,default=[1.,2.])
    p.add_argument("--ring-levels",nargs="+",type=int,default=[2,4,8])
    p.add_argument("--periods",type=int,default=5)
    args = p.parse_args()
    if any(v < 1 for v in args.levels+args.ring_levels) or args.periods < 1:
        p.error("positive levels/periods required")
    if any(v <= 0 or not np.isfinite(v) or abs(2*v-round(2*v))>1e-10 for v in args.lengths):
        p.error("lengths must be positive multiples of .5 m")
    if any(v <= 0 or v > 4 or abs(v-round(v))>1e-10 for v in args.widths):
        p.error("widths must be integer metres in (0,4]")
    out = args.output
    out.mkdir(parents=True,exist_ok=False)
    start = time.perf_counter()
    rings = []
    for ezz in [0.,1e-4]:
        for dim in [2,3]:
            group = []
            for nr in args.ring_levels:
                r = ring(nr,2*nr,dim,ezz)
                if group:
                    factor = nr/group[-1]["nr"]
                    if factor > 1:
                        for q in ["displacement","stress","sigma_zz"]:
                            r[q+"_observed_order"] = float(np.log(group[-1][q+"_relative_L2"]/r[q+"_relative_L2"])/np.log(factor))
                group.append(r)
            rings.extend(group)
        records = []
        for ir,r in enumerate([1.,1.5,2.]):
            for it,t in enumerate([0.,np.pi/4,np.pi/2]):
                ur,sr,st,sz = lame(r,ezz=ezz)
                c,s = np.cos(t),np.sin(t)
                xyz = [r*c,r*s,0.]
                values = [sr*c*c+st*s*s,sr*s*s+st*c*c,sz,(sr-st)*s*c,0.,0.]
                for comp,value in zip(["xx","yy","zz","xy","yz","xz"],values):
                    records.append(record(f"r{ir}-t{it}",xyz,"stress",comp,"Pa",value))
                for comp,value in zip(["x","y","z"],[ur*c,ur*s,0.]):
                    records.append(record(f"r{ir}-t{it}",xyz,"displacement",comp,"m",value))
        name = "lame-plane-strain" if ezz == 0 else "lame-generalized-plane-strain"
        export(out,name,dataset(name,{"E_Pa":30e6,"nu":.25,"a_m":1.,"b_m":2.,
                                      "internal_pressure_Pa":1e5,"external_pressure_Pa":0.,
                                      "axial_strain":ezz,"geometry":"quarter_annulus_symmetry_sides",
                                      "axial_condition":"imposed_uniform_strain_not_closed_end_caps"},records))
    print("Lamé 2D/3D completed",flush=True)
    vib = diagnose(periods=args.periods)
    histories = vib.pop("histories")
    for name,h in histories.items():
        records = []
        for i,t in enumerate(h["time_s"]):
            xyz = [2.,0.,0.] if name == "rod" else [0.,0.,0.]
            for quantity,component,unit,key in [("displacement","x","m","displacement_m"),
                    ("velocity","x","m/s","velocity_m_s"),("acceleration","x","m/s^2","acceleration_m_s2"),
                    ("energy","scalar","J","energy_J")]:
                records.append(record("tip" if name == "rod" else "mass",xyz,quantity,component,unit,h[key][i],t,"monitor"))
        conditions = {"zeta":.03,"initial_tip_displacement_m":.01,"initial_tip_velocity_m_s":.02,
                      "load":"free_response"}
        if name == "rod":
            conditions.update(E_Pa=2e9,rho_kg_m3=2000.,A_m2=.01,L_m=2.,
                              boundary="fixed_free_axial",initial_shape="sin_pi_x_over_2L_tip_normalized",
                              damping="C=2*zeta*continuous_omega1*M")
        else:
            conditions.update(mass_kg=2.,stiffness_N_m=400.,damping="viscous_C=2*zeta*omega*M")
        export(out,name+"-history",dataset(name+"-history",conditions,records,True))
        # Continuous analytical response at exactly the same physical times.
        # Rod temporal error above is against its semidiscrete mode; this
        # second export exposes the continuous PDE spatial+temporal difference.
        analytical_omega = 250*np.pi if name == "rod" else np.sqrt(200.)
        ua,va,aa = exact(np.array(h["time_s"]),analytical_omega)
        analytical_records = []
        modal_mass = 20. if name == "rod" else 2.
        analytical_energy = .5*modal_mass*(va*va+analytical_omega**2*ua*ua)
        for i,t in enumerate(h["time_s"]):
            for quantity,component,unit,values in [("displacement","x","m",ua),
                    ("velocity","x","m/s",va),("acceleration","x","m/s^2",aa),
                    ("energy","scalar","J",analytical_energy)]:
                analytical_records.append(record("tip" if name == "rod" else "mass",xyz,quantity,component,unit,values[i],t,"monitor"))
        analytical_conditions = dict(conditions)
        export(out,name+"-analytic-history",dataset(name+"-history",analytical_conditions,analytical_records,True))
    print("Vibration completed",flush=True)
    bases, boxes = [],[]
    for B in args.widths:
        for level in args.levels:
            pe = foundation(B,level=level,dim=2)
            bases.append(pe)
            for L in args.lengths:
                r = foundation(B,L,level)
                centre = r["profile"][0]["sigma_zz_Pa"]
                limit = pe["profile"][0]["sigma_zz_Pa"]
                r["plane_strain_same_mesh_sigma_zz_Pa"] = limit
                r["central_limit_difference_Pa"] = centre-limit
                r["end_minus_centre_sigma_zz_Pa"] = r["profile"][-1]["sigma_zz_Pa"]-centre
                records = []
                for i,pr in enumerate(r["profile"]):
                    for component in ["zz","yy"]:
                        records.append(record("station-"+str(i),[pr["x_m"],pr["y_m"],pr["z_m"]],
                                              "stress",component,"Pa",pr["sigma_"+component+"_Pa"]))
                name = f"foundation-B{B:g}-L{L:g}-h{level}"
                export(out,name,dataset(name,{"E_Pa":30e6,"nu":.25,"B_m":B,"L_m":L,
                    "domain_x_m":[-2.,2.],"domain_y_m":[-2.,0.],"load_pressure_Pa":1e5,
                    "boundary":"bottom_all_fixed_x_sides_ux_fixed_z_ends_free",
                    "load":"top_uniform_strip_full_length_z",
                    "evaluation":"Q4_Hex8_basis_fixed_points_plus_x_y_trace"},records))
                boxes.append(r)
                print(name+" completed",flush=True)
    mesh = []
    for B in args.widths:
        for L in args.lengths:
            group = sorted([r for r in boxes if r["B_m"]==B and r["L_m"]==L],key=lambda r:r["level"])
            for a,b in zip(group,group[1:]):
                mesh.append({"B_m":B,"L_m":L,"levels":[a["level"],b["level"]],
                             "central_sigma_zz_difference_Pa":b["profile"][0]["sigma_zz_Pa"]-a["profile"][0]["sigma_zz_Pa"],
                             "profile_sigma_zz_RMS_difference_Pa":float(np.sqrt(np.mean([
                                 (bb["sigma_zz_Pa"]-aa["sigma_zz_Pa"])**2 for aa,bb in zip(a["profile"],b["profile"])]))),
                             "order":None,"reason":"two_mesh_levels_no_order_or_acceptance"})
    summary = {"schema":"p6-independent-reference-summary/1","accepted":None,"p6_03_complete":False,
               "original_conditions_reproduced":False,"rings":rings,"vibration":vib,
               "plane_strain_foundation":bases,"finite_foundation":boxes,"mesh_differences":mesh,
               "runtime_seconds":time.perf_counter()-start,"python":platform.python_version(),
               "numpy":np.__version__,"scipy":scipy.__version__,
               "source_sha256":{name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in ["elastic_reference.py","vibration_reference.py","run_reference.py","compare_outputs.py"]}}
    write_json(out/"summary.json",summary)
    print("diagnostic only; accepted=null; P6-03 remains incomplete",flush=True)


if __name__ == "__main__":
    main()
