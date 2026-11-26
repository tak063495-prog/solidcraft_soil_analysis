#!/usr/bin/env bash
# Diagnostic-only run. Archive all solver output, including budget stops.
set -euo pipefail
w=$(realpath "${1:?prepared work}")
d=$(realpath "${2:?output directory}")
mode=${3:?windows, baseline, predicted, or alternate size}
p7=interim/2026-11-26-remaining-p7
p8=interim/2026-11-26-remaining-p8
b="$w/probe/target/release/examples/c01_methods"
export P7_SWITCH_HALF=1
case "$mode" in
    baseline|windows|predicted)
        export P8_RESUME="$w/p9-prefix.json"
        if [[ $mode == windows ]]; then export P9_WINDOWS=1; fi
        if [[ $mode == predicted ]]; then export P8_PREDICT=1; fi
        input="$p7/c01/original1024/input.json"
        ;;
    256|512|1024)
        export P9_HALVES=1
        input="$p7/c01/alternate$mode/input.json"
        ;;
    *) exit 2 ;;
esac
mkdir -p "$d"
python interim/2026-11-25-remaining-p4/profile.py "$b" \
  examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
  "$input" "$d/$mode.json" > "$d/$mode-resource.log" 2> "$w/$mode.trace"
gzip -c "$w/$mode.trace" > "$d/$mode-trace.jsonl.gz"
