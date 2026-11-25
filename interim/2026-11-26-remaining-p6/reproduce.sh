#!/usr/bin/env bash
# Run from repository root in the locked Python/Rust environment.
set -euo pipefail
w=$(realpath -m "${1:?specify a fresh directory}")
mkdir "$w"
p=interim/2026-11-26-remaining-p6
profile=interim/2026-11-25-remaining-p4/profile.py
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir "$w/checks"
bash "$p/prepare.sh" "$w" > "$w/checks/build.log" 2>&1
for level in 1 2; do
    python "$profile" python "$p/face.py" "$w/c02/local$level" "$level" > "$w/checks/local$level.log" 2>&1
done
python "$profile" python "$p/face.py" "$w/c02/focused" 2 --focused > "$w/checks/focused.log" 2>&1
OPENBLAS_NUM_THREADS=8 python "$p/score_face.py" "$w/c02"
for n in 256 512 1024; do
    python "$profile" python "$p/allocation.py" "$w/c01/n$n" "$w/probe/target/release/examples/c01_methods" "$n" > "$w/checks/c01-n$n.log" 2>&1
done
python "$p/score_allocation.py" "$w/c01"
python "$profile" python "$p/nc.py" "$w" "$w/nc" > "$w/checks/nc.log" 2>&1
python "$p/score_nc.py" "$w/nc"
