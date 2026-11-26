#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-26-remaining-p11
base=interim/2026-11-26-remaining-p10
# Insert p11 instrumentation into the archived scratch-only probe build.
sed '/^[(]cd .*cargo build/i python interim/2026-11-26-remaining-p11/instrument.py "$w/probe"' "$base/prepare.sh" > "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
