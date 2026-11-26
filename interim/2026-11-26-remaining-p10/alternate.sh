#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
d=$(realpath "$2")
n=$3
exec 9>"$w/alternate$n.lock"
flock 9
if [[ -f "$w/raw/alternate$n/charges.json" && -f "$w/raw/alternate$n/trace.jsonl.gz" ]]; then exit 0; fi
export P10_KNOTS="$d/plan.json"
python interim/2026-11-25-remaining-p4/profile.py python interim/2026-11-26-remaining-p7/contact.py \
    "$w/probe/target/release/examples/c01_methods" "$w/raw/alternate$n" "$n" alternate > "$d/$n-resource.log"
