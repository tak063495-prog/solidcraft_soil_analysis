#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?empty scratch directory}")
p=interim/2026-11-28-remaining-p15
p14=interim/2026-11-27-remaining-p14
# Preserve the p14 source patch sequence. Compile the new probe only after p15 edits.
sed -e '/^bash "\$w\/prepare-probe.sh"/d' "$p14/prepare.sh" > "$w/generate.sh"
bash "$w/generate.sh" "$w"
sed -i '/^[(]cd .*cargo build/i python interim/2026-11-28-remaining-p15/instrument.py "$w/probe"' "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
