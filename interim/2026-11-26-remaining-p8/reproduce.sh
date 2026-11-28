#!/usr/bin/env bash
# From the repository root, with the locked Python and Rust environment.
set -euo pipefail
w=$(realpath -m "${1:?supply a new work directory}")
mkdir "$w"
out="$w/evidence"
p=interim/2026-11-26-remaining-p8
p4=interim/2026-11-25-remaining-p4
mkdir -p "$out/checks" "$out/c01" "$out/c02" "$out/nc" "$w/solver"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
profile() { python "$p4/profile.py" "$@"; }
# Both reference solves start before either wait. Each uses the original native clock.
pair() {
    profile python "$p/face_pair.py" solve "$out/c02/$1" "$2" > "$out/checks/face-$1.log" 2>&1
    profile python "$p/face_pair.py" residual "$out/c02/$1" > "$out/checks/residual-$1.log" 2>&1
}
pair plus 1 &
pos=$!
pair minus -1 &
neg=$!
wait "$pos"
wait "$neg"
profile python "$p/native_face.py" "$w/solver" "$out/c02" > "$out/checks/native.log" 2>&1
profile python "$p/physical_bands.py" "$w/solver" "$out/c02" > "$out/checks/physical-bands.log" 2>&1
python "$p/face_score.py" "$out/c02/native-active.npz" "$out/c02"
python "$p/alternate.py" "$out/c01/alternate"
bash "$p/contact_runs.sh" "$w/solver" "$out"
python "$p/contact_score.py" "$w/solver" "$out/c01"
profile python "$p/nc.py" "$w/solver" "$out/nc" > "$out/checks/nc.log" 2>&1
