#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?scratch directory}")
p=interim/2026-11-28-remaining-p18
q=interim/2026-11-28-remaining-p17
sed '/^bash "\$w\/prepare-probe.sh"/d' "$q/prepare.sh" > "$w/p17-generate.sh"
bash "$w/p17-generate.sh" "$w"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p18/instrument.py "$w/probe"' "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
