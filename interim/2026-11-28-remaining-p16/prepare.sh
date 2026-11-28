#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?scratch directory}")
p=interim/2026-11-28-remaining-p16
p15=interim/2026-11-28-remaining-p15
sed '/^bash "\$w\/prepare-probe.sh"/d' "$p15/prepare.sh" > "$w/p15-generate.sh"
bash "$w/p15-generate.sh" "$w"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p16/instrument.py "$w/probe"' "$w/prepare-probe.sh"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p16/nc_instrument.py "$w/probe"' "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
