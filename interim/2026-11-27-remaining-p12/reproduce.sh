#!/usr/bin/env bash
# Start in the repository root, in the locked Rust and Python environments.
set -euo pipefail
w=$(realpath -m "${1:?fresh work directory}")
mkdir "$w"
p=interim/2026-11-27-remaining-p12
b=interim/2026-11-26-remaining-p11
export P12_RESULTS="$w/results"
export P11_RESULTS="$P12_RESULTS"
d=$P12_RESULTS
mkdir -p "$d"/{c01/{reload,alternate},c02,nc,checks}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$p/reload.py" freeze "$d/c01/reload"
python "$p/alternate_plan.py" "$d/c01/alternate"
cmp "$d/c01/reload/plan.json" "$p/c01/reload/plan.json"
cmp "$d/c01/alternate/plan.json" "$p/c01/alternate/plan.json"
bash "$b/prepare.sh" "$w" > "$d/checks/build-probe.log" 2>&1
timeout 86400 bash "$p/reload.sh" "$w"
for n in 256 512 1024; do timeout 86400 bash "$b/alternate.sh" "$w" "$n"; done
python "$p/alternate_score12.py" "$w/raw" "$d/c01/alternate"
prof=interim/2026-11-25-remaining-p4/profile.py
for sign in plus minus; do python "$prof" python "$p/face.py" "$sign" "$d/c02" > "$d/checks/face-$sign.log"; done
python "$prof" python "$p/nc.py" "$w" "$d/nc" > "$d/checks/nc.log"
python "$p/verify.py" "$d"
python "$p/resources.py" "$d"
