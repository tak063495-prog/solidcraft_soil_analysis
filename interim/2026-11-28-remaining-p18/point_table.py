"""Render the pointwise reload history and join to the p17 return evidence."""
import json
from pathlib import Path
import sys


def span(pair):
    return "—" if not pair else f"[{pair[0]:.13g}, {pair[1]:.13g}]"


def render(out):
    data=json.loads((out/"summary.json").read_text())
    lines=["# reload から return への点別接線寄与","",
           f"6履歴で一致する reload の端は入口を含め{data['common_endpoints']}。変位・相対変位は m²、反力は kN²。各区間はその工程の進行度。",""]
    for route,title in (("outer","外側"),("middle","中水準")):
        result=data["routes"][route]
        lines.extend([f"## {title}","",
             "| 場・点 | reload入口 | 初増加区間 | 最後の非正→正 | reload終端直前の共通端 | reload終端＝return入口 | return 最初の共通端 |",
             "|---|---:|---|---|---:|---:|---:|"])
        rows=[]
        for k,label in (("u","変位"),("jump","相対変位"),("reaction","反力")):
            for i in data["points"]:
                point=result[k][str(i)]
                s,j=point["tangent"],point["join"]
                rows.append((label,i,s,j))
                lines.append(f"| {label}・{i} | {s['at_start']:+.7e} | {span(s['first_increase'])} | "
                             f"{span(j['last_nonpositive_then_positive'])} | {j['penultimate']:+.7e} | "
                             f"{j['at_reload_end']:+.7e} | {j['at_return_first']:+.7e} |")
        lines.extend(["","| 場・点 | reload の最大増加区間 | 増加量 | return の最大増加区間 | 増加量 |",
                      "|---|---|---:|---|---:|"])
        for label,i,s,j in rows:
            lines.append(f"| {label}・{i} | {span(s['largest_increase'])} | {s['max_increase']:+.7e} | "
                         f"{span(j['return_max_span'])} | {j['return_max_increase']:+.7e} |")
        native=rows[0][-1]["reload_last_half"]
        lines.extend(["","終端に到着した各水準の実際の受理半区間は "+
                      "、".join(f"N{n} {span(t)}" for n,t in native.items())+"。",""])
    (out/"points.ja.md").write_text("\n".join(lines).rstrip()+"\n")


if __name__=="__main__":
    render(Path(sys.argv[1]))
