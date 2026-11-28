"""Audit subsequent material substeps on the two sides of each actual partition."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"2026-11-26-remaining-p9"))
from nc_scan import P8,p8


def run(work,out):
    p=p8.prior.old
    probe=work/"probe/target/release/examples"
    m=p.nc_extended.Extended(work/"inputs",probe/"p5_homogeneous")
    old=p.read(P8/"nc/parent-632.json.gz")
    states,lam,t=old["active_material_states"],old["load"],np.array(old["secant"])
    main=json.loads((out/"summary.json").read_text())
    targets={s["parent"]:s for s in main["parents"]}
    ordinary=[r for r in map(json.loads,(out/"global.jsonl").read_text().splitlines()) if r.get("control") in ("replay","ordinary")]
    checked=set()
    results=[]
    calls=global_calls=0
    with gzip.open(out/"substeps.jsonl.gz","xt") as f:
        for row in ordinary:
            index=row["parent"]
            assert lam==row["load"]
            if index in targets and index not in checked:
                checked.add(index)
                flat=np.zeros_like(m.base)
                flat[m.free]=t[:-1]
                slope=-np.einsum("eqai,ei->eqa",m.b,flat[m.dofs]).reshape(-1,6)
                for i in (272,351):
                    cut=m.width/np.abs(slope[i]).sum()
                    assert cut==targets[index]["partition_boundaries"][str(i)]
                    cache={}
                    def sample(mult):
                        nonlocal calls
                        if mult not in cache:
                            strain=(cut*mult)*slope[i]
                            needed=max(1,int(np.ceil(np.abs(strain).sum()/m.width)))
                            n=2*needed if needed>1 else 1
                            rec=dict(point=i,law=m.laws[m.ids[i]],state=states[i],strain=strain.tolist(),count=n)
                            try:
                                obs=p8.point(probe,rec)
                            except AssertionError as error:
                                # Preserve P8's cross-check. An ill-conditioned sample
                                # cannot certify a sign or shrink a bracket.
                                obs=dict(steps=[],result=dict(reason="P8 Schur/direction audit: "+str(error)))
                            calls+=1
                            # A substep may not be reached if an earlier return fails.
                            ds={s["step"]:s["denominator"] for s in obs["steps"]}
                            cache[mult]=dict(mult=float(mult),width=cut*mult,count=n,denominators=ds)
                            print(json.dumps(dict(parent=index,point=i,**cache[mult],disposition=obs["result"]["reason"]),
                                             allow_nan=False),file=f)
                        return cache[mult]
                    grid=sorted(set(np.r_[np.linspace(.05,.99,40),.995,.999,.99999,1-1e-9,1+1e-9,1.00001,1.001,1.005,np.linspace(1.01,1.99,40)]))
                    vals=[sample(x) for x in grid]
                    brackets=[]
                    for a0,b0 in zip(vals,vals[1:]):
                        if a0["count"]!=b0["count"]:
                            continue
                        for j in a0["denominators"].keys() & b0["denominators"].keys():
                            a,b=a0,b0
                            if a["denominators"][j]*b["denominators"][j]>=0:
                                continue
                            interrupted=False
                            for _ in range(14):
                                mid=sample((a["mult"]+b["mult"])/2)
                                if j not in mid["denominators"]:
                                    interrupted=True
                                    break
                                if a["denominators"][j]*mid["denominators"][j]>0:
                                    a=mid
                                else:
                                    b=mid
                            brackets.append(dict(substep=j,left=a,right=b,interrupted=interrupted))
                    results.append(dict(parent=index,load=lam,point=i,cut=cut,left=sample(1-1e-9),
                                        right=sample(1+1e-9),same_partition_brackets=brackets))
            if index==max(targets):
                break
            res,it,why=p.nc_branch.correct(m,states,lam,t,row["width"],1e-12)
            global_calls+=1
            assert (res is not None)==row["accepted"] and why==row["reason"]
            if res is not None:
                dz,states,_=res
                lam+=dz[-1]
                t=dz/np.linalg.norm(dz)
    failed_log=out/"substeps-initial.jsonl.gz"
    failed_global=failed_local=0
    if failed_log.exists():
        with gzip.open(failed_log,"rt") as src:
            prior=[json.loads(l) for l in src]
        failed_local=len(prior)+1
        last=prior[-1]["parent"]
        failed_global=sum(r["parent"]<last for r in ordinary)
    total=main["total_including_p9"]+failed_global+global_calls
    assert len(checked)==6 and total<=2000
    (out/"substeps.json").write_text(json.dumps(dict(parents=results,extra_global=global_calls,extra_local=calls,
                                                      initial_run_global=failed_global,initial_run_local=failed_local,total_global_including_p9=total),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
