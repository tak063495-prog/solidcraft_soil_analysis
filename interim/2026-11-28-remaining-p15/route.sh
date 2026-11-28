#!/usr/bin/env bash
set -euo pipefail
p=interim/2026-11-28-remaining-p15
p14=interim/2026-11-27-remaining-p14
p7=interim/2026-11-26-remaining-p7
w=$(realpath "$1")
d=${P15_RESULTS:-$p}
prof=interim/2026-11-25-remaining-p4/profile.py
bin="$w/probe/target/release/examples/c01_methods"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
if [[ "$2" == reload ]]; then
    dest="$d/c01/reload"
    mkdir -p "$dest"
    P8_PREFIX="$w/prefix.json" P7_SWITCH_HALF=1 timeout 86400 python "$prof" "$bin" \
       examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
       "$p7/c01/original1024/input.json" "$w/prefix-result.json" \
       > "$dest/prefix-resource.log" 2> "$w/prefix.trace"
    gzip -nc "$w/prefix.trace" > "$dest/prefix.jsonl.gz"
    python interim/2026-11-26-remaining-p9/prefix_compare.py \
       interim/2026-11-26-remaining-p10/c01/reload/reverse-end.json.gz "$w/prefix.json" \
       "$w/prefix.trace" "$dest/prefix-identity.json"
    P8_RESUME="$w/prefix.json" P7_SWITCH_HALF=1 P10_LABELS=1 P14_GAP=1 P15_RELOAD=1 \
       P11_CROSS=interim/2026-11-27-remaining-p12/c01/reload/plan.json \
       timeout 86400 python "$prof" "$bin" examples/joint-water-column-tet10-mesh.json \
         examples/ground-joint-water-mc.json "$p7/c01/original1024/input.json" "$dest/result.json" \
         > "$dest/resource.log" 2> "$w/reload.trace"
    gzip -nc "$w/reload.trace" > "$dest/trace.jsonl.gz"
    python "$p/reload_score.py" "$dest"
else
    n="$2"
    dest="$d/c01/alternate/outer"
    mkdir -p "$dest"
    P11_OBSERVE=1 P15_OUTER=1 P10_KNOTS="$p14/c01/alternate/outer/plan.json" \
      timeout 86400 python "$prof" python "$p7/contact.py" "$bin" "$w/outer/alternate$n" "$n" alternate \
         > "$dest/$n-resource.log" 2>&1
fi
