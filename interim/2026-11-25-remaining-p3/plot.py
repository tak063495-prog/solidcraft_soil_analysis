"""Plot audited histories. Usage: python plot.py EVIDENCE_FOLDER OUTPUT.png."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot(root, output):
    c02=json.loads((root/"c02/summary.json").read_bytes())["body_to_fixed_reference"]
    fig,(a,b)=plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for key,label,style in [("slip_relative","Area-weighted slip","o-"),
                            ("signed_area_relative","Signed area integral","s--")]:
        a.plot([3,4,5],[100*c02[f"p{p}"][key] for p in (3,4,5)],style,label=label)
    a.axhline(0.1,color="0.35",linestyle=":",label="Adopted 0.1% width")
    a.set(xticks=[3,4,5],xticklabels=["P3 / 256","P4 / 512","P5 / 1024"],yscale="log",
          xlabel="Displacement degree / base increments",ylabel="Relative L2 difference (%)",
          title="C-02 · Fixed independent reference")
    a.legend(frameon=False,fontsize=9)
    folders=["0.32","0.04","0.02","unit-0.04","roundtrip-0.04"]
    for i,folder in enumerate(folders):
        path=root/"nc"/folder
        meta=json.loads((path/"summary.json").read_bytes())
        rows=[json.loads(s) for s in (path/"trials.jsonl").read_text().splitlines()]
        rows=[r for r in rows if r["accepted"]]
        arc=[0.]+[r["arc_m"] for r in rows]
        load=[meta["start_lambda"]]+[r["load"] for r in rows]
        label=f"Max step {folder} m" if i<3 else ["Unit check","Round-trip unit check"][i-3]
        line,=b.plot(arc,load,["-","--",":","-.",":"][i],lw=1.5,label=label)
        b.scatter(arc[-1],load[-1],s=18,color=line.get_color())
    b.set(xlabel="Accumulated continuation step (m)",ylabel="Load factor λ",
          title="NC-01 · From the simple fold")
    b.legend(frameon=False,fontsize=8)
    for ax in (a,b):
        ax.grid(alpha=0.2)
    fig.savefig(output,dpi=180)
    plt.close(fig)


if __name__=="__main__":
    plot(Path(sys.argv[1]),Path(sys.argv[2]))
