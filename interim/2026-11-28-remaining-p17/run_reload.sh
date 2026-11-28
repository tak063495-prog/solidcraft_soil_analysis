#!/usr/bin/env bash
# Rebuild and verify the reverse endpoint before running the new reload.
set -euo pipefail
p=interim/2026-11-28-remaining-p16
export P17_WORK=1
bash "$p/reload.sh" "$1" "$2"
