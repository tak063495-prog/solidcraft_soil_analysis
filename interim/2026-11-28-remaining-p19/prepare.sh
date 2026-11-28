#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?scratch directory}")
p=interim/2026-11-28-remaining-p18
sed '/^bash "\$w\/prepare-probe.sh"/d' "$p/prepare.sh" > "$w/p18-generate.sh"
bash "$w/p18-generate.sh" "$w"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p19/instrument.py "$w/probe"' "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
