#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?work}")
out=$(realpath "${2:?alternate evidence directory}")
n=${3:?level}
p7=interim/2026-11-26-remaining-p7
export P7_SWITCH_HALF=1 P9_MATCHED="$out/matched-plan.json" P9_LEVEL="$n"
python interim/2026-11-25-remaining-p4/profile.py "$w/matched/target/release/examples/c01_methods" \
  examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
  "$p7/c01/alternate$n/input.json" "$out/matched-$n.json" \
  > "$out/matched-$n-resource.log" 2> "$w/matched-$n.trace"
gzip -c "$w/matched-$n.trace" > "$out/matched-$n-trace.jsonl.gz"
