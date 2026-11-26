#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-27-remaining-p12
d=${P12_RESULTS:-$p}/c01/reload
bin="$w/probe/target/release/examples/c01_methods"
prof=interim/2026-11-25-remaining-p4/profile.py
input=interim/2026-11-26-remaining-p7/c01/original1024/input.json
args=(examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json "$input")
export P7_SWITCH_HALF=1
python "$prof" env P8_PREFIX="$w/prefix.json" "$bin" "${args[@]}" "$w/prefix-result.json" > "$d/prefix-resource.log" 2> "$w/prefix.trace"
gzip -nc "$w/prefix.trace" > "$d/prefix.jsonl.gz"
python "$p/reload.py" prefix "$d" "$w"
python "$prof" env P8_RESUME="$w/prefix.json" P11_CROSS="$d/plan.json" "$bin" "${args[@]}" "$d/result.json" > "$d/entrance-resource.log" 2> "$w/entrance.trace"
gzip -nc "$w/entrance.trace" > "$d/entrance.jsonl.gz"
python "$p/reload.py" score "$d"
