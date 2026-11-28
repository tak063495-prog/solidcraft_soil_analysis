#!/usr/bin/env bash
# From the repository root, in the locked Rust / Python environment.
set -euo pipefail
w=$(realpath -m "${1:?supply a fresh output directory}")
mkdir "$w"
p=interim/2026-11-26-remaining-p7
p4=interim/2026-11-25-remaining-p4
p6=interim/2026-11-26-remaining-p6
mkdir "$w/checks"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
profile() { python "$p4/profile.py" "$@"; }
profile python "$p/face_shift.py" "$w/c02" > "$w/checks/c02.log" 2>&1
OPENBLAS_NUM_THREADS=8 python "$p/face_shift.py" "$w/c02" --score
profile python "$p6/face.py" "$w/c02/repeat" 2 --focused > "$w/checks/c02-repeat.log" 2>&1
python "$p/face_shift.py" "$w/c02" --repeat "$w/c02/repeat"
bash "$p6/prepare.sh" "$w" > "$w/checks/build-nc.log" 2>&1
profile python "$p/nc_pair.py" "$w" "$w/nc" > "$w/checks/nc.log" 2>&1
# Replace the P6 contact-only instrumentation with the P7 opt-in experiment.
cp src/stages_joint_path_error.rs "$w/probe/src/"
python "$p/prepare.py" "$w/probe"
(cd "$w/probe" && cargo build --locked --release --example c01_methods) > "$w/checks/build-contact.log" 2>&1
for route in original alternate; do
    for n in 256 512 1024; do
        profile python "$p/contact.py" "$w/probe/target/release/examples/c01_methods" \
            "$w/raw/$route$n" "$n" "$route" > "$w/checks/$route$n.log" 2>&1
    done
done
python "$p/contact.py" score "$w/raw" "$w/c01"
