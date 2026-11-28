"""Replay the original small B-02/B-03 cases and their unchanged comparators."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
B02 = ROOT / "interim/2026-12-01-tri6-offset-composite-b02"
B03 = ROOT / "interim/2026-12-02-tri6-offset-kink-b03"
sys.path[:0] = [str(B02), str(B03)]
from cases import inputs, read, write
from b03_cases import make
from compare import compare
from study import command, exact, sha
from scope import run as run_b02_scope
from probe import run as run_b03_scope
from supplemental import run as run_b03_supplement


def runs_for(label):
    if label == "B-02":
        for use_beam in (False, True):
            for sign,e in (("negative",-0.05),("zero",0.),("positive",0.05)):
                name = ("beam-" if use_beam else "plate-") + sign
                yield name,inputs(0,4,e,use_beam),4
        for steps in (8,16):
            yield f"load-{steps}",inputs(0,steps,0.05,True),steps
    else:
        for kind in ("planar","kink-zero","normal-sign-pair","step-split","kink-split","height-step"):
            for steps in ((4,) if kind=="planar" else (4,8,16)):
                yield f"{kind}-h0-s{steps}",make(kind,0,steps),steps


def evaluate(suite,output,binary,cli):
    results = []
    source = B02 if suite=="B-02" else B03
    history = {r["case"]:r for r in read(source/"measurements.json")["cases"]}
    for name,(mesh,ground,stages,_),steps in runs_for(suite):
        # One original path at a time. Keep comparisons, inputs and process
        # logs; the full field streams are regenerable from the listed hashes.
        case_out=output/name
        case_out.mkdir(parents=True)
        with tempfile.TemporaryDirectory(dir=output,delete=False) as temp:
            folder=Path(temp)
            base=[]
            for key,value in (("mesh",mesh),("ground",ground),("stages",stages)):
                file=folder/(key+".json")
                write(file,value)
                # Match the original Linux input blob on Windows as well.
                file.write_bytes(file.read_bytes().replace(b"\r\n",b"\n"))
                (case_out/file.name).write_bytes(file.read_bytes())
                base.append(file)
            hashes={p.name:sha(p) for p in base}
            if hashes != history[name]["input_sha256"]:
                raise ValueError(f"archived input hash mismatch for {name}")
            stops = [2*steps+steps//2]
            if name=="kink-zero-h0-s4":
                stops=[steps//2,steps+steps//2,*stops,3*steps+steps//2,4*steps+steps//2]
            if name in ("beam-negative","beam-positive"):
                stops+=[steps+steps//2,3*steps+steps//2,4*steps+steps//2]
            native=command([binary,"solve",*base,folder/"native",",".join(map(str,stops))],case_out/"native.log")
            command([sys.executable,B02/"reference.py",folder,folder/"independent.jsonl"],case_out/"independent.log")
            comparison=compare(folder)
            (case_out/"comparison.json").write_bytes((folder/"comparison.json").read_bytes())
            command([binary,"plain",*base,folder/"plain"],case_out/"plain.log")
            original=read(folder/"native/run.json")
            if not original["complete"] or not exact(original,read(folder/"plain/run.json")):
                raise ValueError(f"observed/unobserved mismatch for {name}")
            accepted=[json.loads(s) for s in (folder/"native/accepted.jsonl").read_text().splitlines()]
            resumed_count=0
            for stop in stops:
                target=folder/f"resumed-{stop}"
                cp=folder/f"native/checkpoint-{stop}.json"
                resumed=command([binary,"resume",*base,target,cp],case_out/f"resume-{stop}.log")
                later=[json.loads(s) for s in (target/"accepted.jsonl").read_text().splitlines()]
                if native["pid"]==resumed["pid"] or not exact(original,read(target/"run.json")) or not exact(later,[s for s in accepted if s["accepted_index"]>stop]):
                    raise ValueError(f"checkpoint mismatch {name} {stop}")
                resumed_count+=len(later)
            attempts=original["geostatic"]["attempts"]+original["attempts"]
            accepted_attempts=[item for item in attempts if item["accepted"]]
            force=max(item["relative_force_residual"] for item in accepted_attempts)
            linear=max((entry["true_relative_residual"] for item in accepted_attempts for entry in item["linear_solves"]),default=0.)
            with (folder/"independent.jsonl").open() as stream:
                direct=max(json.loads(line).get("direct_linear_residual",0.) for line in stream)
            gates=[comparison["existing_field_gates_passed"],comparison["maximum_native_reassembled_residual"]<1e-10,
                   comparison["maximum_independent_reassembled_residual"]<1e-10]
            # Preserve the inclusive B-02 and strict B-03 comparisons in the original collectors.
            if suite=="B-02":
                gates.extend([force<=1e-8,linear<=1e-10,direct<1e-10])
            else:
                gates.extend([force<1e-8,linear<1e-10])
            results.append(dict(case=name,inputs_sha256=hashes,tet10=len(mesh["tetrahedra"]),
                                bounded_gates_passed=all(gates),force_residual=force,true_linear_residual=linear,
                                independent_direct_linear_residual=direct,
                                accepted_points=comparison["accepted_points"],checkpoint_indices=stops,
                                later_states_compared=resumed_count,existing_field_gates_passed=comparison["existing_field_gates_passed"],
                                native_residual=comparison["maximum_native_reassembled_residual"],
                                independent_residual=comparison["maximum_independent_reassembled_residual"]))
        shutil.rmtree(temp)
        write(output/"fields.json",results)
        print(name,results[-1]["accepted_points"],flush=True)
    if suite=="B-02":
        run_b02_scope(output/"scope",cli)
        probes=read(output/"scope/refusals.json")["records"]
    else:
        probes=run_b03_scope(output/"scope",cli)["records"]
        run_b03_supplement(output/"supplemental",cli)
    compact=[]
    for row in probes:
        compact.append({k:row[k] for k in ("name","category","refused","refused_as_expected","accepted","complete","matched_expected","diagnostic","incomplete_scope_response") if k in row})
    write(output/"results.json",dict(suite=suite,passed=False,fields=results,scope=compact))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("suite",choices=["B-02","B-03"])
    parser.add_argument("output",type=Path)
    filename=lambda name:name+(".exe" if os.name=="nt" else "")
    parser.add_argument("--observer",type=Path,default=ROOT/"target/release"/filename("b02-observer"))
    parser.add_argument("--cli",type=Path,default=ROOT/"target/release"/filename("soil-fem-3d"))
    args=parser.parse_args()
    evaluate(args.suite,args.output.resolve(),args.observer.resolve(),args.cli.resolve())
