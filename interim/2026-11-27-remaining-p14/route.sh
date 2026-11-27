#!/usr/bin/env bash
# From the repository root after prepare.sh, with a fresh scratch directory.
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-27-remaining-p14
prof=interim/2026-11-25-remaining-p4/profile.py
bin="$w/probe/target/release/examples/c01_methods"
d=${P14_RESULTS:-$p}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ "$2" == reload ]]; then
    dest="$d/c01/reload"
    P8_PREFIX="$w/prefix.json" P7_SWITCH_HALF=1 timeout 86400 python "$prof" "$bin" \
        examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
        interim/2026-11-26-remaining-p7/c01/original1024/input.json "$w/prefix-result.json" \
        > "$dest/prefix-resource.log" 2> "$w/prefix.trace"
    gzip -nc "$w/prefix.trace" > "$dest/prefix.jsonl.gz"
    python interim/2026-11-26-remaining-p9/prefix_compare.py \
        interim/2026-11-26-remaining-p10/c01/reload/reverse-end.json.gz "$w/prefix.json" \
        "$w/prefix.trace" "$dest/prefix-identity.json"
    P8_RESUME="$w/prefix.json" P7_SWITCH_HALF=1 P10_LABELS=1 P14_GAP=1 \
        P11_CROSS=interim/2026-11-27-remaining-p12/c01/reload/plan.json \
        timeout 86400 python "$prof" "$bin" examples/joint-water-column-tet10-mesh.json \
        examples/ground-joint-water-mc.json interim/2026-11-26-remaining-p7/c01/original1024/input.json \
        "$dest/result.json" > "$dest/resource.log" 2> "$w/reload.trace"
    gzip -nc "$w/reload.trace" > "$dest/trace.jsonl.gz"
else
    route="$2"
    n="$3"
    dest="$d/c01/alternate/$route"
    P11_OBSERVE=1 P10_KNOTS="$dest/plan.json" timeout 86400 python "$prof" \
        python interim/2026-11-26-remaining-p7/contact.py "$bin" "$w/$route/alternate$n" \
        "$n" alternate > "$dest/$n-resource.log" 2>&1
fi
