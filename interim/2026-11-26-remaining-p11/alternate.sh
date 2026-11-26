#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
n=$2
p=interim/2026-11-26-remaining-p11
d="${P11_RESULTS:-$p}/c01/alternate"
export P11_OBSERVE=1 P10_KNOTS="$d/plan.json"
python interim/2026-11-25-remaining-p4/profile.py python interim/2026-11-26-remaining-p7/contact.py \
    "$w/probe/target/release/examples/c01_methods" "$w/raw/alternate$n" "$n" alternate > "$d/$n-resource.log"
