#!/usr/bin/env bash
# Run from the repository root with the project's pinned Rust and Python environment.
set -euo pipefail
p=interim/2026-11-28-remaining-p18
p16=interim/2026-11-28-remaining-p16
prof=interim/2026-11-25-remaining-p4/profile.py
w=$(realpath -m "${1:?fresh scratch directory}")
mkdir "$w"
out="$w/results"
mkdir -p "$out"/{checks,nc,c02,c01/{reload/control,alternate}}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
bash "$p/prepare.sh" "$w" > "$out/checks/build.log" 2>&1
for sign in minus plus; do
    timeout 86400 python "$prof" python "$p/face.py" "$sign" "$out/c02" > "$out/checks/face-$sign.log" 2>&1
done
mkdir "$w/face-independent"
for sign in minus plus; do
    timeout 86400 python "$prof" python "$p/face.py" "$sign" "$w/face-independent" > "$out/checks/face-independent-$sign.log" 2>&1
done
python interim/2026-11-28-remaining-p17/check_face_replay.py "$out/c02" "$w/face-independent" "$out/checks/face-replay.json"
mkdir "$out/c02/minus-328"
timeout 86400 python "$prof" python "$p/face.py" minus "$out/c02/minus-328" 328 > "$out/checks/face-minus-328.log" 2>&1
timeout 86400 bash "$p/original.sh" "$w" "$out/c01/reload"
timeout 86400 bash "$p/original.sh" "$w" "$out/c01/reload/control" control
python "$p/score_reload.py" "$out/c01/reload"
for route in outer middle; do
    for n in 256 512 1024; do
        timeout 86400 python "$prof" python "$p/reload_halves.py" "$w" "$out/c01/alternate" "$route" "$n" > "$out/checks/$route-$n.log" 2>&1
    done
done
python "$p16/check_replays.py" "$w" "$out/c01/alternate"
python "$p/reload_split.py" "$w" "$out/c01/alternate"
python "$p/point_table.py" "$out/c01/alternate"
python "$p16/gates.py" "$out/c01/alternate"
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
timeout 86400 python "$prof" python "$p/nc_bracket.py" "$w" "$out/nc" > "$out/checks/nc.log" 2>&1
P18_RESULTS="$out" python "$p/verify.py" > "$out/checks/p18-tests.log" 2>&1
cargo test --release --locked --test joint_law --test joint_path_error --test mcc_midpoint > "$out/checks/rust.log" 2>&1
for n in {9..17}; do
    scripts=(interim/*remaining-p"$n"/verify.py)
    python "${scripts[0]}" > "$out/checks/p$n-tests.log" 2>&1
done
git rev-parse 21b2530 > "$out/checks/base.txt"
python "$p/resources.py" "$w" "$out"
