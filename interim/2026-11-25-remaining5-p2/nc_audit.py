"""Reassemble every saved equilibrium using the existing full-state audit."""
import gzip
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
from nc_continue import HERE, Section, prepare, internal


def audit(inputs,folder):
    section=Section(inputs,None)
    mesh=json.loads((inputs/"mesh.json").read_bytes())
    with zipfile.ZipFile(HERE.parent/"2026-11-23-nc01-coordinate-polish/evidence.zip") as archive:
        run=json.loads(archive.read("run.json"))
    cp=run["last_committed"]
    q,off=prepare(mesh,run["integration_points"])
    target_u=np.ravel(cp["displacement_m"])
    fixed=np.setdiff1d(np.arange(len(target_u)),section.free)
    output={}
    for path in sorted(folder.glob("0.*")):
        report=json.loads((path/"report.json").read_bytes())
        assert report["start"]==section.start and not report["stop"]
        data=report["rows"]
        assert len(data)==report["requested_steps"]==len(report["attempts"])
        ratios=[]
        yields=[]
        with gzip.open(path/"accepted.jsonl.gz","rt") as f:
            for row,line in zip(data,f,strict=True):
                state=json.loads(line)
                assert (state["step"],state["load"])==(row["step"],row["load"])
                np.testing.assert_array_equal(np.ravel(state["displacement_m"])[fixed],target_u[fixed])
                materials=list(cp["material_states"])
                for i,point in zip(section.point_ids,state["active_material_states"],strict=True):
                    materials[i]=point
                force,y=internal(mesh,run["ground"],materials,cp["active_elements"],q,off,
                                 material_ids=cp["element_material_ids"])
                load=section.load(row["load"])
                ratios.append(float(np.linalg.norm((load-force.ravel())[section.free])/max(1.,np.linalg.norm(load))))
                yields.append(float(y))
        assert max(ratios)<=section.settings["relative_force_tolerance"]
        linear=[it for a in report["attempts"] for it in a["iterations"] if "augmented_true_linear" in it]
        lr=max(max(i["augmented_true_linear"],i["force_true_linear"]) for i in linear)
        assert lr<=section.settings["relative_linear_tolerance"]
        index=int(np.argmax([r["load"] for r in data]))
        assert 0<index<len(data)-1 and all(r["delta_lambda"]>0 for r in data[:index+1])
        assert all(r["delta_lambda"]<0 for r in data[index+1:])
        output[path.name]=dict(accepted=len(data),positive=index+1,negative=len(data)-index-1,
            peak_lambda=data[index]["load"],peak_arc_m=data[index]["arc_m"],end_lambda=data[-1]["load"],
            max_independent_relative_force=max(ratios),max_normalized_yield=max(yields),
            max_original_linear_rows_or_augmented=lr)
    return output


if __name__=="__main__":
    print(json.dumps(audit(Path(sys.argv[1]),Path(sys.argv[2])),indent=2,allow_nan=False))
