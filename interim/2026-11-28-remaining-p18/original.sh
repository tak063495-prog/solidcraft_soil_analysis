#!/usr/bin/env bash
# Resume from the identical reverse endpoint used in p17. This is a solver checkpoint.
set -euo pipefail
w=$(realpath "$1")
dest=$(realpath "$2")
mode=${3:-shadow}
run="$w/reload-$mode"
mkdir -p "$run"
unset P18_CONTROL
case "$mode" in
    shadow) ;;
    control) export P18_CONTROL=1 ;;
    *) exit 2 ;;
esac
p10=interim/2026-11-26-remaining-p10
p12=interim/2026-11-27-remaining-p12
p7=interim/2026-11-26-remaining-p7
common="$dest"
if [[ "$mode" == control ]]; then common=$(dirname "$dest"); fi
(
    flock -x 9
    if [[ ! -f "$common/prefix-identity.json" ]]; then
        P8_PREFIX="$w/prefix.json" P7_SWITCH_HALF=1 timeout 86400 python interim/2026-11-25-remaining-p4/profile.py \
          "$w/probe/target/release/examples/c01_methods" examples/joint-water-column-tet10-mesh.json \
          examples/ground-joint-water-mc.json "$p7/c01/original1024/input.json" "$w/prefix-result.json" \
          > "$common/prefix-resource.log" 2> "$w/prefix.trace"
        gzip -nc "$w/prefix.trace" > "$common/prefix.jsonl.gz"
        python interim/2026-11-26-remaining-p9/prefix_compare.py "$p10/c01/reload/reverse-end.json.gz" \
          "$w/prefix.json" "$w/prefix.trace" "$common/prefix-identity.json"
    fi
) 9>"$w/prefix.lock"
P8_RESUME="$w/prefix.json" P7_SWITCH_HALF=1 P10_LABELS=1 P14_GAP=1 P15_RELOAD=1 \
 P11_CROSS="$p12/c01/reload/plan.json" P16_SLIP=1 P17_WORK=1 P18_OBSERVE=1 \
 timeout 86400 python interim/2026-11-25-remaining-p4/profile.py \
 "$w/probe/target/release/examples/c01_methods" examples/joint-water-column-tet10-mesh.json \
 examples/ground-joint-water-mc.json "$p7/c01/original1024/input.json" "$dest/result.json" \
 > "$dest/resource.log" 2> "$run/original.trace"
gzip -nc "$run/original.trace" > "$dest/trace.jsonl.gz"
