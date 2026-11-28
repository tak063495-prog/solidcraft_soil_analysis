#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?scratch directory}")
p=interim/2026-11-28-remaining-p17
p16=interim/2026-11-28-remaining-p16
sed '/^bash "\$w\/prepare-probe.sh"/d' "$p16/prepare.sh" > "$w/p16-generate.sh"
bash "$w/p16-generate.sh" "$w"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p17/instrument.py "$w/probe"' "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
