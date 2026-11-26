#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
n=$2
p=interim/2026-11-26-remaining-p11
d="${P11_RESULTS:-$p}/c01/alternate"
p10=interim/2026-11-26-remaining-p10
export P7_SWITCH_HALF=1 P9_HALVES=1 P11_OBSERVE=1 P10_KNOTS="$p10/c01/alternate/plan.json"
b="$w/probe/target/release/examples/c01_methods"
python interim/2026-11-25-remaining-p4/profile.py "$b" examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
 "$p10/c01/alternate/$n/input.json" "$w/observed-$n.json" > "$d/observe-$n-resource.log" 2> "$w/observed-$n.trace"
gzip -nc "$w/observed-$n.trace" > "$d/observed-$n.jsonl.gz"
