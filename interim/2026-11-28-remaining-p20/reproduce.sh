#!/usr/bin/env bash
# From the repository root with its pinned toolchain and numerical packages.
set -euo pipefail
w=$(realpath -m "${1:?fresh work directory}")
mkdir "$w"
out="$w/results"
mkdir -p "$out"/{checks,c01/{original,alternate},c02,nc}
p=interim/2026-11-28-remaining-p20
prof=interim/2026-11-25-remaining-p4/profile.py
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
bash "$p/prepare.sh" "$w" > "$out/checks/build.log" 2>&1
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original > "$out/checks/inputs.log"
for sign in minus plus; do
    python "$prof" python "$p/face_refine.py" "$sign" "$out/c02" > "$out/checks/face-$sign.log" 2>&1
done
python "$p/summarize_face.py" "$out/c02"
python "$prof" python "$p/alternate.py" "$out/c01/alternate" > "$out/checks/alternate.log" 2>&1
python "$prof" python interim/2026-11-28-remaining-p16/gates.py "$out/c01/alternate" > "$out/checks/alternate-gates.log" 2>&1
python "$prof" python "$p/original.py" solve "$w" "$out/c01/original" > "$out/checks/original.log" 2>&1
gzip -nc "$w/original.trace" > "$out/c01/original/trace.jsonl.gz"
python "$p/original.py" score "$w" "$out/c01/original"
python "$p/original.py" slip "$w" "$out/c01/original"
python "$p/slip_summary.py" "$out/c01/original"
for n in 256 512; do
    python "$prof" python interim/2026-11-26-remaining-p7/contact.py "$w/probe/target/release/examples/c01_methods" "$w/recheck/original$n" "$n" original > "$out/checks/original$n.log" 2>&1
done
python "$p/recheck.py" "$w" "$out/c01/original"
python "$prof" python "$p/nc.py" "$w" "$out/nc" > "$out/checks/nc.log" 2>&1
python "$p/nc_summary.py" "$out/nc"
P20_RESULTS="$out" python "$p/verify.py"
python interim/2026-11-28-remaining-p19/verify.py
cargo test --release --locked --test joint_law --test joint_path_error --test mcc_midpoint
