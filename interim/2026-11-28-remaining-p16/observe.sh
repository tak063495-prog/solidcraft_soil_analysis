#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
out=$(realpath "$2")
route=$3
n=$4
p=interim/2026-11-28-remaining-p16
timeout 86400 python interim/2026-11-25-remaining-p4/profile.py python "$p/observe.py" "$w" "$out/c01/alternate" "$route" "$n" > "$out/checks/$route-$n.log" 2>&1
