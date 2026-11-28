#!/usr/bin/env bash
# Run at the repository root in the pinned Rust/Python environment.
set -euo pipefail
p=interim/2026-11-28-remaining-p19
prof=interim/2026-11-25-remaining-p4/profile.py
w=$(realpath -m "${1:?fresh scratch directory}")
mkdir "$w"
out="$w/results"
mkdir -p "$out"/{checks,c01/{reload/control,alternate},nc,c02/{baseline,tight}}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
bash "$p/prepare.sh" "$w" > "$out/checks/build.log" 2>&1
for variant in baseline tight; do
    signs=(minus plus)
    opts=()
    if [[ "$variant" == tight ]]; then signs=(minus); opts=(--tight); fi
    for s in "${signs[@]}"; do
        timeout 86400 python "$prof" python "$p/face.py" "$s" "$out/c02/$variant" "${opts[@]}" \
            > "$out/checks/$variant-$s.log" 2>&1
    done
done
python "$p/leaf_score.py" "$out/c02"
bash "$p/reload.sh" "$w" "$out/c01/reload" prefix
bash "$p/reload.sh" "$w" "$out/c01/reload" shadow
bash "$p/reload.sh" "$w" "$out/c01/reload/control" control
python "$p/point19.py" "$out/c01/reload"
python "$p/score_reload.py" "$out/c01/reload"
timeout 86400 python "$prof" python "$p/alternate.py" "$out/c01/alternate" > "$out/checks/alternate.log" 2>&1
timeout 86400 python "$prof" python interim/2026-11-28-remaining-p16/gates.py "$out/c01/alternate" > "$out/checks/gates.log" 2>&1
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
timeout 86400 python "$prof" python "$p/nc_next.py" "$w" "$out/nc" > "$out/checks/nc.log" 2>&1
P19_RESULTS="$out" python "$p/verify.py"
for n in {9..18}; do
    scripts=(interim/*remaining-p"$n"/verify.py)
    python "${scripts[0]}"
done
cargo test --release --locked --test joint_law --test joint_path_error --test mcc_midpoint
