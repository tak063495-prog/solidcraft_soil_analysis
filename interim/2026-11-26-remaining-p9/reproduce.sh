#!/usr/bin/env bash
# Run from repository root in the locked Rust/Python environment.
set -euo pipefail
w=$(realpath -m "${1:?supply a fresh work directory}")
mkdir "$w"
d="$w/evidence"
mkdir -p "$d"/{c01/{alternate,reload},c02,nc,checks}
p=interim/2026-11-26-remaining-p9
p8=interim/2026-11-26-remaining-p8
f=interim/2026-11-25-remaining-p4/profile.py
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
python "$p/freeze_windows.py" > "$w/windows.json"
cmp "$w/windows.json" "$p/windows.json"
bash "$p/prepare.sh" "$w" > "$d/checks/build.log" 2>&1
for sign in plus minus; do
    python "$f" python "$p/face_bounds.py" "$sign" "$d/c02" > "$d/checks/bounds-$sign.log" 2>&1
done
python "$f" python "$p/equilibrium.py" "$d/c02" > "$d/checks/equilibrium.log" 2>&1
python "$p/score_face.py" "$d/c02"
P7_SWITCH_HALF=1 P8_PREFIX="$w/p9-prefix.json" python "$f" "$w/probe/target/release/examples/c01_methods" \
    examples/joint-water-column-tet10-mesh.json examples/ground-joint-water-mc.json \
    interim/2026-11-26-remaining-p7/c01/original1024/input.json "$d/c01/reload/prefix.json" \
    > "$d/checks/prefix-resource.log" 2> "$w/prefix.trace"
python "$p/prefix_compare.py" "$p8/c01/reverse-end-checkpoint.json.gz" "$w/p9-prefix.json" \
    "$w/prefix.trace" "$d/c01/reload/prefix-comparison.json"
for mode in windows baseline predicted; do
    bash "$p/contact_run.sh" "$w" "$d/c01/reload" "$mode"
done
python "$p/score_reload.py" "$d/c01/reload"
for n in 256 512 1024; do
    bash "$p/contact_run.sh" "$w" "$d/c01/alternate" "$n"
done
python "$p/alternate.py" "$d/c01/alternate"
python "$p/within_halves.py" "$d/c01/alternate"
cp -a "$w/probe" "$w/matched"
python "$p/matched_prepare.py" "$d/c01/alternate" "$w/matched/src"
(cd "$w/matched" && cargo build --release --locked --example c01_methods) > "$d/checks/matched-build.log" 2>&1
for n in 256 512 1024; do
    bash "$p/matched_run.sh" "$w" "$d/c01/alternate" "$n"
done
python "$p/matched_score.py" "$d/c01/alternate"
python "$f" python "$p/nc_scan.py" "$w" "$d/nc" > "$d/checks/nc.log" 2>&1
python "$f" python "$p/nc_partition.py" "$w" "$d/nc" > "$d/checks/nc-partition.log" 2>&1

python "$p/resources.py" "$d"
