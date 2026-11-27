#!/usr/bin/env bash
# From repository root, with the project's locked Python and Rust tools on PATH.
set -euo pipefail
p=interim/2026-11-28-remaining-p15
w=$(realpath -m "${1:?new scratch directory}")
mkdir "$w"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export P15_RESULTS="$w/results"
d=$P15_RESULTS
mkdir -p "$d"/{checks,c02,nc,c01/{reload,alternate/outer}}
bash "$p/prepare.sh" "$w" > "$d/checks/build.log" 2>&1
python "$p/schur.py" "$d/nc/schur.json"
for sign in plus minus; do
    timeout 86400 python interim/2026-11-25-remaining-p4/profile.py python "$p/face_direct.py" "$sign" "$d/c02" > "$d/checks/face-$sign.log" 2>&1
done
bash "$p/route.sh" "$w" reload
for n in 256 512 1024; do
    bash "$p/route.sh" "$w" "$n"
done
python "$p/alternate_score_p15.py" "$w" "$d"
python "$p/verify.py" > "$d/checks/p15-tests.log" 2>&1
python "$p/resources.py" "$d"
cargo test --locked --release --test joint_law --test joint_path_error --test mcc_midpoint > "$d/checks/rust.log" 2>&1
for n in 9 10 11 12 13 14; do
    script=(interim/*remaining-p"$n"/verify.py)
    python "${script[0]}" > "$d/checks/p$n-tests.log" 2>&1
done
