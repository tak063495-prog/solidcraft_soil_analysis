#!/usr/bin/env bash
set -euo pipefail
p15=interim/2026-11-28-remaining-p15
w=$(realpath "$1")
export P15_RESULTS=$(realpath "$2")
export P16_SLIP=1
mkdir -p "$P15_RESULTS/c01/reload"
sed '/^    python "\$p\/reload_score.py"/d' "$p15/route.sh" > "$w/p16-route.sh"
timeout 86400 bash "$w/p16-route.sh" "$w" reload
