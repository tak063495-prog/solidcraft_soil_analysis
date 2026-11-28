"""Human-readable stage and point comparisons from the additive arrays."""
import json
from pathlib import Path
import sys

import numpy as np

STAGES = ("geostatic","shear","open","recontact","reverse","reload","return")
KEYS = ("plastic_slip_dot","accumulated_dot","u_point_growth","jump_point_growth","reaction_point_growth")
LABEL = ("塑性滑り内積","累積滑り内積","変位二乗差","相対変位二乗差","反力二乗差")
UNITS = ("m²","m²","m²","m²","kN²")


def table(out):
    lines = ["# 接触点と工程の照合","",
        "内積は隣り合う水準差の積、二乗差は細分後の差の二乗ノルムから細分前を引いた値。工程表の二乗差は実測値で、点別表は接触トラクションによる寄与である。各荷重・各場の未説明分も NPZ に保存した。",""]
    stages = json.loads((out/"attribution.json").read_text())
    for route,title in (("outer","外側"),("middle","中水準")):
        lines += ["## "+title,"","| 工程 | "+" | ".join(s+" ("+u+")" for s,u in zip(LABEL,UNITS))+" |",
                  "|---|"+"---:|"*len(KEYS)]
        for stage in STAGES:
            row = stages[route]["stages"][stage]
            numbers = [sum(row[k]) for k in KEYS[:2]]
            for k in ("u","jump","reaction"):
                coarse,fine = row[k+"_squared_norms"]
                numbers.append(fine-coarse)
            lines.append("| "+stage+" | "+" | ".join(f"{v:.5e}" for v in numbers)+" |")
        with np.load(out/f"{route}-1792.npz") as f:
            for item,sl in (("再接触1.0",slice(4*256-1,4*256)),("return 全256荷重",slice(6*256,None))):
                values = [f[k][sl].sum(0) for k in KEYS]
                chosen = set(np.argsort(values[0])[:3])
                for v in values[2:]:
                    chosen.update(np.argsort(-abs(v))[:3])
                order = sorted(chosen,key=lambda i:values[0][i])
                lines += ["","### "+item,"","塑性滑り内積が最も負の3点と、各力学量の寄与の絶対値が最大の3点の和集合。","",
                    "| 接触点 | "+" | ".join(s+" ("+u+")" for s,u in zip(LABEL,UNITS))+" |",
                    "|---:|"+"---:|"*len(KEYS)]
                for i in order:
                    lines.append(f"| {i} | "+" | ".join(f"{v[i]:.5e}" for v in values)+" |")
                lines.append("| 32点計 | "+" | ".join(f"{v.sum():.5e}" for v in values)+" |")
        lines.append("")
    (out/"points-and-stages.ja.md").write_text("\n".join(lines).rstrip()+"\n")


if __name__=="__main__":
    table(Path(sys.argv[1]))
