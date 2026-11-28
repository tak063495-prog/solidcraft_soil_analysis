#!/usr/bin/env bash
# In the locked environment, run in the repository root with a fresh work directory.
set -euo pipefail
p=interim/2026-11-26-remaining-p5
profile=interim/2026-11-25-remaining-p4/profile.py
w=$(realpath -m "${1:?specify a new directory}")
mkdir "$w"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
bash "$p/prepare.sh" "$w"
for d in 4 8; do
    python "$profile" python "$p/c02_reference.py" "$w/c02/d$d" "$d"
done
python "$p/c02_compare.py" "$w/c02"
for n in 256 512 1024; do
    python "$profile" python "$p/c01_trace.py" "$w/c01/n$n" "$w/probe/target/release/examples/c01_methods" "$n"
done
python "$p/c01_analyze.py" "$w/c01"
python "$profile" python "$p/nc_trace.py" "$w" "$w/nc"
python "$p/nc_analyze.py" "$w/nc"
cargo test --locked --release --test joint_law --test joint_path_error --test mcc_midpoint
