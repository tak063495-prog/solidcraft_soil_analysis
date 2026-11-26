#!/usr/bin/env bash
# Run from the repository root with the locked toolchain activated.
set -euo pipefail
w=$(realpath "${1:?work directory}")
p=interim/2026-11-26-remaining-p8
d=$(realpath "${2:?evidence directory}")
bash interim/2026-11-26-remaining-p6/prepare.sh "$w" > "$d/checks/build.log" 2>&1
cp src/stages_joint_path_error.rs "$w/probe/src/"
python interim/2026-11-26-remaining-p7/prepare.py "$w/probe"
python "$p/prepare_contact.py" "$w/probe"
(cd "$w/probe" && cargo build --locked --release --example c01_methods) >> "$d/checks/build.log" 2>&1
mkdir "$w/parents"
export P7_SWITCH_HALF=1
b="$w/probe/target/release/examples/c01_methods"
profile=interim/2026-11-25-remaining-p4/profile.py
args=(examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json interim/2026-11-26-remaining-p7/c01/original1024/input.json)
P8_PREFIX="$w/prefix.json" python "$profile" "$b" "${args[@]}" "$d/c01/prefix.json" > "$w/prefix.out" 2> "$w/prefix.trace"
export P8_RESUME="$w/prefix.json"
P8_PREDICT=1 python "$profile" "$b" "${args[@]}" "$d/c01/predicted.json" > "$w/predicted.out" 2> "$w/predicted.trace" &
pred_job=$!
P8_SAVE_PARENTS="$w/parents" python "$profile" "$b" "${args[@]}" "$d/c01/baseline.json" > "$w/baseline.out" 2> "$w/baseline.trace"
P8_PREDICT=1 P8_REPLAY="$w/parents" python "$profile" "$b" "${args[@]}" "$d/c01/paired.json" > "$w/paired.out" 2> "$w/paired.trace"
wait "$pred_job"
