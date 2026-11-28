#!/usr/bin/env bash
# Run at repository root with the locked Python and Rust environments.
set -euo pipefail
work=$(realpath -m "${1:?new work directory}")
mkdir "$work"
p=interim/2026-11-27-remaining-p13
export P13_RESULTS="$work/results"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
d=$P13_RESULTS
mkdir -p "$d"/{c01/{alternate,reload},c02,nc,checks}
python "$p/plans.py" "$d/c01/alternate"
bash interim/2026-11-26-remaining-p11/prepare.sh "$work" > "$d/checks/build.log" 2>&1
prof=interim/2026-11-25-remaining-p4/profile.py
for sign in plus minus; do
    timeout 86400 python "$prof" python "$p/local_face.py" "$sign" "$d/c02" > "$d/checks/face-$sign.log" 2>&1
done
timeout 86400 bash "$p/reload.sh" "$work"
for point in 14 26; do
    for n in 256 512 1024; do
        bash "$p/alternate.sh" "$work" "$point" "$n"
    done
done
python "$p/alternate.py" "$work" "$d/c01/alternate"
timeout 86400 python "$prof" python "$p/continue_nc.py" "$work" "$d/nc" > "$d/checks/nc.log" 2>&1
python "$p/nc_check.py" "$work" "$d/nc"
python "$p/verify.py" "$d" > "$d/checks/p13-tests.log" 2>&1
python "$p/resources.py" "$d"
