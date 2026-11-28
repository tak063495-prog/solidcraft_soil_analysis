#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
d=$(realpath "$2")
p9=interim/2026-11-26-remaining-p9
p7=interim/2026-11-26-remaining-p7
p=$(dirname "$(realpath "$0")")
export P7_SWITCH_HALF=1
b="$w/probe/target/release/examples/c01_methods"
run() {
    local label=$1
    shift
    python interim/2026-11-25-remaining-p4/profile.py env "$@" "$b" \
        examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
        "$p7/c01/original1024/input.json" "$d/$label.json" > "$d/$label-resource.log" 2> "$w/$label.trace"
    gzip -c "$w/$label.trace" > "$d/$label-trace.jsonl.gz"
}
run prefix P8_PREFIX="$w/prefix.json"
gzip -nc "$w/prefix.json" > "$d/reverse-end.json.gz"
python "$p/prefix_identity.py" "$w/prefix.json" "$w/prefix.trace" "$d/prefix-identity.json"
export P8_RESUME="$w/prefix.json"
run observation P9_WINDOWS=1 P10_LABELS=1
python "$p/freeze.py" "$d/observation-trace.jsonl.gz" > "$d/plan.json"
sha256sum "$d/plan.json" | cut -d' ' -f1 > "$d/plan.sha256"
run narrow P10_WINDOWS="$d/plan.json"
