#!/usr/bin/env bash
# Run at the repository root with the locked Python and Rust environments.
set -euo pipefail
w=$(realpath -m "${1:?new scratch directory}")
mkdir "$w"
p=interim/2026-11-26-remaining-p10
d="$w/evidence"
mkdir -p "$d"/{c01/{alternate,reload},c02,nc,checks}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
prof=interim/2026-11-25-remaining-p4/profile.py
python "$p/freeze.py" > "$d/c01/alternate/plan.json"
cmp "$p/c01/alternate/plan.json" "$d/c01/alternate/plan.json"
bash "$p/prepare.sh" "$w" > "$d/checks/build.log" 2>&1
bash "$p/reload.sh" "$w" "$d/c01/reload"
python "$p/score_contact.py" reload "$d/c01/reload"
for n in 256 512 1024; do
    bash "$p/alternate.sh" "$w" "$d/c01/alternate" "$n"
done
python "$p/score_contact.py" alternate "$w/raw" "$d/c01/alternate"
for s in plus minus; do
    python "$prof" python "$p/face.py" "$s" "$d/c02" > "$d/checks/face-$s-resource.log"
    python "$prof" python "$p/tangent_samples.py" "$s" "$d/c02" > "$d/checks/tangent-$s-resource.log"
done
python "$p/score_face.py" "$d/c02"
python "$prof" python "$p/nc.py" "$w" "$d/nc" > "$d/checks/nc-resource.log"
python "$prof" python "$p/nc_substeps.py" "$w" "$d/nc" > "$d/checks/nc-substeps-resource.log"
python "$p/verify.py" "$d"
python "$p/resources.py" "$d"
