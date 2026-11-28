#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-26-remaining-p11
d="${P11_RESULTS:-$p}/c01/reload"
p10=interim/2026-11-26-remaining-p10
export P7_SWITCH_HALF=1
b="$w/probe/target/release/examples/c01_methods"
input=interim/2026-11-26-remaining-p7/c01/original1024/input.json
python interim/2026-11-25-remaining-p4/profile.py env P8_PREFIX="$w/prefix.json" "$b" \
    examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json "$input" "$w/reload-prefix-summary.json" \
    > "$d/prefix-resource.log" 2> "$w/reload-prefix.trace"
gzip -nc "$w/reload-prefix.trace" > "$d/prefix.jsonl.gz"
python "$p/prefix_check.py" "$w/prefix.json" "$w/reload-prefix.trace" "$d/prefix-identity.json"
export P8_RESUME="$w/prefix.json"
for mode in cross ordinary half; do
    unset P10_WINDOWS P11_CROSS
    case "$mode" in
        cross) export P11_CROSS="$d/plan.json" ;;
        half) export P10_WINDOWS="$p10/c01/reload/plan.json" ;;
    esac
    python interim/2026-11-25-remaining-p4/profile.py "$b" examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
        interim/2026-11-26-remaining-p7/c01/original1024/input.json "$w/reload-$mode.json" \
        > "$d/$mode-resource.log" 2> "$w/reload-$mode.trace"
    gzip -nc "$w/reload-$mode.trace" > "$d/$mode.jsonl.gz"
done
