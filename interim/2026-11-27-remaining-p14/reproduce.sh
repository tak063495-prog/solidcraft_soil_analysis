#!/usr/bin/env bash
# From the requested repository root, with the locked Rust/Python tools on PATH.
set -euo pipefail
p=interim/2026-11-27-remaining-p14
w=$(realpath -m "${1:?new empty work directory}")
mkdir "$w"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export P14_RESULTS="$w/results"
d=$P14_RESULTS
mkdir -p "$d"/{checks,c02,c01/{alternate,reload},nc}
python "$p/plans.py" "$d/c01/alternate"
bash "$p/prepare.sh" "$w" > "$d/checks/build.log" 2>&1
for sign in plus minus; do
    timeout 86400 python interim/2026-11-25-remaining-p4/profile.py python "$p/face_chain.py" "$sign" "$d/c02" > "$d/checks/face-$sign.log"
done
python interim/2026-11-25-remaining-p4/profile.py python "$p/local_lu.py" "$w/probe/target/release/examples/local_lu" "$d/nc" > "$d/checks/local-lu.log"
timeout 86400 bash "$p/route.sh" "$w" reload
python "$p/reload_score.py" "$d/c01/reload"
for route in middle outer; do
    for n in 256 512 1024; do
        bash "$p/route.sh" "$w" "$route" "$n"
    done
done
python "$p/alternate.py" "$w" "$d/c01/alternate"
python "$p/verify.py" "$d" > "$d/checks/p14-tests.log" 2>&1
python "$p/resources.py" "$d"
