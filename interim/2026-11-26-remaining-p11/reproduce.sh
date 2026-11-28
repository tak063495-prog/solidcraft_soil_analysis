#!/usr/bin/env bash
# Run from the repository root in the locked Rust / Python environment.
set -euo pipefail
w=$(realpath -m "${1:?fresh scratch directory}")
mkdir "$w"
p=interim/2026-11-26-remaining-p11
p10=interim/2026-11-26-remaining-p10
export P11_RESULTS="$w/evidence"
d="$P11_RESULTS"
mkdir -p "$d"/{c01/{alternate,reload},c02/weighted,nc,checks}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
prof=interim/2026-11-25-remaining-p4/profile.py
python "$p/crossings.py" "$p10/c01/reload/narrow-trace.jsonl.gz" "$p10/c01/reload/plan.json" > "$d/c01/reload/plan.json"
cmp "$p/c01/reload/plan.json" "$d/c01/reload/plan.json"
bash "$p/prepare.sh" "$w" > "$d/checks/prepare.log" 2>&1
bash "$p/reload.sh" "$w"
python "$p/reload_score.py" "$d/c01/reload"
for n in 256 512 1024; do bash "$p/observe.sh" "$w" "$n"; done
python "$p/alt_plan.py" "$d/c01/alternate"
cmp "$p/c01/alternate/plan.json" "$d/c01/alternate/plan.json"
for n in 256 512 1024; do bash "$p/alternate.sh" "$w" "$n"; done
python "$p/alternate_score.py" "$w/raw" "$d/c01/alternate"
python "$p/matched.py" "$d/c01/alternate"
for sign in plus minus; do
    for radius in 0 1e-8; do
        python "$prof" env P11_TEST_WEIGHT=1 python "$p/face_matrix.py" "$sign" "$d/c02" "$radius" > "$d/checks/face-$sign-$radius.log"
        python "$prof" env P11_TEST_WEIGHT=4 python "$p/face_matrix.py" "$sign" "$d/c02/weighted" "$radius" > "$d/checks/face-weighted-$sign-$radius.log"
    done
done
python "$p/face_score.py" "$d/c02"
python "$prof" python "$p/nc_paths.py" "$w" "$d/nc" > "$d/checks/nc.log"
python "$p/verify.py" "$d"
python "$p/resources.py" "$d"
