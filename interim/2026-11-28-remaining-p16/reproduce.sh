#!/usr/bin/env bash
# Run from the repository root with the locked Python and Rust toolchains.
set -euo pipefail
w=$(realpath -m "${1:?new scratch directory}")
mkdir "$w"
p=interim/2026-11-28-remaining-p16
out="$w/results"
mkdir -p "$out"/{checks,nc,c02,c01/{alternate,reload}}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
prof=interim/2026-11-25-remaining-p4/profile.py
bash "$p/prepare.sh" "$w" > "$out/checks/build.log" 2>&1
for sign in plus minus; do
    timeout 86400 python "$prof" python "$p/paired_face.py" "$sign" "$out/c02" > "$out/checks/face-$sign.log" 2>&1
done
python "$p/reload_cost.py" "$out/c01/reload"
bash "$p/reload.sh" "$w" "$out"
python "$p/reload_cost.py" "$out/c01/reload"
for route in middle outer; do
    for n in 256 512 1024; do bash "$p/observe.sh" "$w" "$out" "$route" "$n"; done
done
python "$p/check_replays.py" "$w" "$out/c01/alternate"
python "$p/attribution.py" "$out/c01/alternate" --dense
python "$p/gates.py" "$out/c01/alternate"
python "$p/field_tables.py" "$out/c01/alternate"
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
timeout 86400 python "$prof" python "$p/parent_nc.py" "$w" "$out/nc" > "$out/checks/nc.log" 2>&1
python "$p/resources.py" "$w" "$out"
P16_RESULTS="$out" python "$p/verify.py"
cargo test --locked --release --test joint_law --test joint_path_error --test mcc_midpoint > "$out/checks/rust.log" 2>&1
for n in 9 10 11 12 13 14 15; do
    script=(interim/*remaining-p"$n"/verify.py)
    python "${script[0]}" > "$out/checks/p$n-tests.log" 2>&1
done
