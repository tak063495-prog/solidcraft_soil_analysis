#!/usr/bin/env bash
set -euo pipefail
work=$(realpath "$1")
point=$2
level=$3
p=interim/2026-11-27-remaining-p13
dest=${P13_RESULTS:-$p}/c01/alternate/point$point
export P11_OBSERVE=1 P10_KNOTS="$dest/plan.json"
timeout 86400 python interim/2026-11-25-remaining-p4/profile.py \
    python interim/2026-11-26-remaining-p7/contact.py \
    "$work/probe/target/release/examples/c01_methods" "$work/point$point/alternate$level" \
    "$level" alternate > "$dest/$level-resource.log" 2>&1
