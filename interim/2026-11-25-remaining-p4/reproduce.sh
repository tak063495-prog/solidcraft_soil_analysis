#!/usr/bin/env bash
# Run in the repository root in the locked Rust/Python environment.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
p=interim/2026-11-25-remaining-p4
p3=interim/2026-11-25-remaining-p3
w=${1:?specify a new working directory}
mkdir "$w"
w=$(realpath "$w")
cargo build --locked --release --lib --bin soil-fem-3d --example c01_methods --message-format=json > "$w/artifacts.jsonl"
serde=$(python - "$w/artifacts.jsonl" <<'PY'
import json,sys
print(next(x['filenames'][0] for x in map(json.loads,open(sys.argv[1])) if x.get('target',{}).get('name')=='serde_json'))
PY
)
args=(--edition 2024 -O -Dwarnings -L dependency=target/release/deps --extern soil_fem_3d=target/release/libsoil_fem_3d.rlib --extern "serde_json=$serde")
for spec in '4 256 1' '5 256 1' '5 512 2' '3 512 2'; do
    read -r d n split <<< "$spec"
    tag=p${d}n${n}
    python "$p/profile.py" python "$p3/c02_run.py" "$w/raw-$tag" "$w/c02/$tag" target/release/soil-fem-3d "$d" "$n" "$split"
done
python "$p/c02_diagnose.py" "$w/c02"
for spec in '256 original' '512 original' '1024 original' '256 alternate'; do
    read -r n route <<< "$spec"
    python "$p/profile.py" python "$p/c01_probe.py" "$w/$route$n" "$n" "$route"
done
python "$p/c01_diagnose.py" "$w" "$w/c01"
python "$p/profile.py" python "$p/c01_local.py" > "$w/c01-local.jsonl"
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
rustc "${args[@]}" interim/2026-11-25-remaining5-p2/nc_map.rs -o "$w/nc_map"
sed 's/powi(448)/powi(480)/' "$p3/nc_units.rs" > "$w/nc480.rs"
rustc "${args[@]}" --cfg roundtrip_units "$p3/nc_units.rs" -o "$w/nc448"
rustc "${args[@]}" --cfg roundtrip_units "$w/nc480.rs" -o "$w/nc480"
rustc "${args[@]}" "$p/nc_homogeneous.rs" -o "$w/nc_homogeneous"
rustc "${args[@]}" "$p/nc_point.rs" -o "$w/nc_point"
mkdir "$w/nc"
for spec in 'fine nc_map 0.02' 'aux nc448 roundtrip-0.04'; do
    read -r name worker source <<< "$spec"
    python "$p/profile.py" python "$p/nc_secant.py" "$w/inputs" "$w/$worker" "$w/replay-$name" "$source"
    cmp "$w/replay-$name/secant.npy" "$p/nc/$name-secant.npy"
    cp "$w/replay-$name/secant.npy" "$w/nc/$name-secant.npy"
done
for spec in 'fine448 0.02 448' 'aux448 roundtrip-0.04 448' 'aux480 roundtrip-0.04 480'; do
    read -r name source exponent <<< "$spec"
    export NC_WORKER="$w/nc$exponent" NC_FAILURES="$w/nc/$name-incoming-failures.jsonl"
    python "$p/profile.py" python "$p/nc_restart.py" "$w/inputs" "$p/nc_capture.py" "$w/nc/$name-incoming" "$source"
done
"$w/nc_point" < "$w/nc/aux448-incoming-failures.jsonl" > "$w/nc/aux-point-replays.jsonl"
python "$p/nc_verify.py" "$w/inputs" "$w"/nc/*-incoming > "$w/nc/physical-checks.json"
python "$p/profile.py" python "$p/nc_extended.py" "$w/inputs" "$w/nc_homogeneous" "$w/nc/fine-extended-retry" 0.02
python "$p/nc_extended_audit.py" "$w/inputs" "$w/nc448" "$w/nc_homogeneous" "$w/nc/fine-extended-retry" "$w/nc/fine448-incoming" > "$w/nc/extended-audit.json"
cargo test --locked --release --test joint_law --test joint_path_error --test mcc_midpoint
python "$p/profile.py" python "$p/nc_extended.py" "$w/inputs" "$w/nc_homogeneous" "$w/nc/extended-004" 0.02 14200 0.04
python "$p/nc_extended_audit.py" "$w/inputs" "$w/nc448" "$w/nc_homogeneous" "$w/nc/extended-004" "$w/nc/fine448-incoming" > "$w/nc/extended-004-audit.json"
python "$p/nc_compare.py" "$w/nc" > "$w/nc/path-comparison.json"
