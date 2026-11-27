#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-27-remaining-p14
# Apply the archived instrumentation in original order and compile only once.
sed -e '/^python interim\/2026-11-22-impl-linear-solver\/prepare_inputs.py/d' \
    -e 's/--example p5_homogeneous)/--example p5_homogeneous --example local_lu)/' \
    -e '/^[(]cd .*cargo build/i python interim/2026-11-26-remaining-p11/instrument.py "$w/probe"\npython interim/2026-11-27-remaining-p14/instrument.py "$w/probe"\ncp interim/2026-11-27-remaining-p14/local_lu.rs "$w/probe/examples/local_lu.rs"' \
    interim/2026-11-26-remaining-p10/prepare.sh > "$w/prepare-probe.sh"
bash "$w/prepare-probe.sh" "$w"
