#!/usr/bin/env bash
# Run from the repository root in the pinned project environment. Use a fresh scratch path.
set -euo pipefail
p=interim/2026-11-28-remaining-p17
p16=interim/2026-11-28-remaining-p16
w=$(realpath -m "${1:?new scratch directory}")
mkdir "$w"
out="$w/results"
mkdir -p "$out"/{checks,nc,c02,c01/{reload,alternate}}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
prof=interim/2026-11-25-remaining-p4/profile.py
bash "$p/prepare.sh" "$w" > "$out/checks/build.log" 2>&1
for sign in plus minus; do
    python "$prof" python "$p/leaf_face.py" "$sign" "$out/c02" > "$out/checks/face-$sign.log" 2>&1
done
mkdir "$w/face-replay"
for sign in minus plus; do
    python "$prof" python "$p/leaf_face.py" "$sign" "$w/face-replay" > "$out/checks/face-$sign-replay.log" 2>&1
done
python "$p/check_face_replay.py" "$out/c02" "$w/face-replay" "$out/checks/face-replay.json"
python "$p/c02/leaf_totals.py" "$out/c02"
timeout 86400 bash "$p/run_reload.sh" "$w" "$out" > "$out/checks/reload.log" 2>&1
python "$p/reload_score17.py" "$out/c01/reload"
for route in outer middle; do
    for n in 256 512 1024; do
        python "$prof" python "$p/return_history.py" "$w" "$w" "$route" "$n" > "$out/checks/$route-$n.log" 2>&1
    done
done
python "$p16/check_replays.py" "$w" "$out/c01/alternate"
python "$p/return_split.py" "$w" "$out/c01/alternate"
python "$p/return_table.py" "$out/c01/alternate"
python "$p16/gates.py" "$out/c01/alternate"
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
python "$prof" python "$p/nc_parent.py" "$w" "$out/nc" > "$out/checks/nc.log" 2>&1
python "$p/nc_conditions.py" "$out/nc"
P17_RESULTS="$out" python "$p/verify.py" > "$out/checks/p17-tests.log" 2>&1
cargo test --release --locked --test joint_law --test joint_path_error --test mcc_midpoint > "$out/checks/rust.log" 2>&1
for n in 9 10 11 12 13 14 15 16; do
    scripts=(interim/*remaining-p"$n"/verify.py)
    python "${scripts[0]}" > "$out/checks/p$n-tests.log" 2>&1
done
python "$p/resources.py" "$w" "$out"
