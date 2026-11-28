#!/usr/bin/env bash
set -euo pipefail
work=$(realpath "$1")
p=interim/2026-11-27-remaining-p13
dest=${P13_RESULTS:-$p}/c01/reload
python "$p/reload.py" freeze "$dest"
bin="$work/probe/target/release/examples/c01_methods"
prof=interim/2026-11-25-remaining-p4/profile.py
args=(examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json
      interim/2026-11-26-remaining-p7/c01/original1024/input.json)
export P7_SWITCH_HALF=1
timeout 86400 python "$prof" env P8_PREFIX="$work/prefix.json" "$bin" "${args[@]}" "$work/prefix-result.json" \
    > "$dest/prefix-resource.log" 2> "$work/prefix.trace"
gzip -nc "$work/prefix.trace" > "$dest/prefix.jsonl.gz"
python interim/2026-11-27-remaining-p12/reload.py prefix "$dest" "$work"
timeout 86400 python "$prof" env P8_RESUME="$work/prefix.json" P10_LABELS=1 P11_CROSS="$dest/plan.json" \
    "$bin" "${args[@]}" "$dest/result.json" > "$dest/resource.log" 2> "$work/reload.trace"
gzip -nc "$work/reload.trace" > "$dest/trace.jsonl.gz"
python "$p/reload.py" score "$dest"
