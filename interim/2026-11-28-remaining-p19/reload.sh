#!/usr/bin/env bash
# Run from the repository root after prepare.sh. The prefix is common to both trials.
set -euo pipefail
w=$(realpath "$1")
out=$(realpath "$2")
mode=$3
p=interim/2026-11-26-remaining-p7
prof=interim/2026-11-25-remaining-p4/profile.py
bin="$w/probe/target/release/examples/c01_methods"
args=(examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json "$p/c01/original1024/input.json")
export P7_SWITCH_HALF=1
unset P19_CONTROL P18_CONTROL P8_RESUME P8_PREFIX
case "$mode" in
    prefix)
        P8_PREFIX="$w/prefix.json" timeout 86400 python "$prof" "$bin" "${args[@]}" "$w/prefix-result.json" \
            > "$out/prefix-resource.log" 2> "$w/prefix.trace"
        gzip -nc "$w/prefix.trace" > "$out/prefix.jsonl.gz"
        python interim/2026-11-26-remaining-p9/prefix_compare.py \
            interim/2026-11-26-remaining-p10/c01/reload/reverse-end.json.gz \
            "$w/prefix.json" "$w/prefix.trace" "$out/prefix-identity.json"
        exit ;;
    control) export P19_CONTROL=1 ;;
    shadow) ;;
    *) exit 2 ;;
esac
export P8_RESUME="$w/prefix.json" P10_LABELS=1 P14_GAP=1 P15_RELOAD=1
export P11_CROSS=interim/2026-11-27-remaining-p12/c01/reload/plan.json
export P16_SLIP=1 P17_WORK=1 P18_OBSERVE=1 P18_CONTROL=1
timeout 86400 python "$prof" "$bin" "${args[@]}" "$out/result.json" > "$out/resource.log" 2> "$w/$mode.trace"
gzip -nc "$w/$mode.trace" > "$out/trace.jsonl.gz"
