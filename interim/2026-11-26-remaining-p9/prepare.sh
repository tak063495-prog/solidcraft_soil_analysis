#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?new scratch directory}")
p=interim/2026-11-26-remaining-p9
bash interim/2026-11-26-remaining-p6/prepare.sh "$w"
cp src/stages_joint_path_error.rs "$w/probe/src/"
python interim/2026-11-26-remaining-p7/prepare.py "$w/probe"
python interim/2026-11-26-remaining-p8/prepare_contact.py "$w/probe"
python "$p/instrument.py" "$w/probe"
(cd "$w/probe" && cargo build --release --locked --example c01_methods)
