"""Render stage totals and the signed onset from the measured common return endpoints."""
from pathlib import Path
import json
import sys


def interval(span):
    return "なし" if span is None else " → ".join(f"{v:.10g}" for v in span)


def run(out):
    data=json.loads((out/"summary.json").read_text())
    lines=["# return の法線・接線寄与と最初の増加","","点は従来の0始まり。二乗ノルムの隣接水準差に対する寄与は符号付き。工程計は従来の256格子端の和。「初増加」は直前の共通端から値が初めて増えた区間、「初正」は終端値が初めて正になる区間。入口ですでに正の寄与がある場合も入口値を明示した。","","各共通端の直前の受理半区間は NPZ の half_256、half_512、half_1024 に収録した。異なる水準で細かい固有の端が入る場合、表の両端は共通端での識別範囲である。"]
    for route,title in (("outer","外側"),("middle","中水準")):
        lines+=["","## "+title,"","| 場 | 点 | 成分 | 入口値 | 工程計 | 初増加 | 初正 | 最大増加区間 |",
                "|---|---:|---|---:|---:|---|---|---|"]
        for field,label in (("u","変位 m²"),("jump","相対変位 m²"),("reaction","反力 kN²")):
            for point in data["points"]:
                for component,name in (("normal","法線"),("tangent","接線"),("sum","計")):
                    v=data["routes"][route][field][str(point)][component]
                    lines.append(f"| {label} | {point} | {name} | {v['at_start']:.5e} | {v['grid256_sum']:.5e} | "
                        f"{interval(v['first_increase'])} | {interval(v['first_positive'])} | {interval(v['largest_increase'])} |")
    (out/"points.ja.md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":
    run(Path(sys.argv[1]))
