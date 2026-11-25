#!/usr/bin/env bash
# From the repository root, in the locked Python/Rust environment.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
p=interim/2026-11-25-remaining-p3
p2=interim/2026-11-25-remaining5-p2
w=${1:?supply a new output directory}
mkdir "$w"
w=$(realpath "$w")
cargo build --release --locked --lib --bin soil-fem-3d --example c01_methods \
    --message-format=json > "$w/build-artifacts.jsonl"
serde=$(python - "$w/build-artifacts.jsonl" <<'PY'
import json,sys
for line in open(sys.argv[1]):
    item=json.loads(line)
    if item.get("target",{}).get("name")=="serde_json":
        print(next(f for f in item["filenames"] if f.endswith(".rlib")))
        break
else:
    raise SystemExit("serde_json rlib not found")
PY
)
for spec in '3 256 1' '4 512 2' '5 1024 2'; do
    read -r degree n sub <<< "$spec"
    python "$p/c02_run.py" "$w/c02-raw/p$degree" "$w/c02/p$degree" \
        target/release/soil-fem-3d "$degree" "$n" "$sub"
done
python "$p/c02_score.py" "$w/c02"
for n in 256 512 1024; do
    python "$p/c01_route.py" "$w/fixed/n$n" "$n" --fixed-budget
done
for n in 256 512; do
    python "$p/c01_route.py" "$w/scaled/n$n" "$n"
done
ln -s "$w/fixed/n1024" "$w/scaled/n1024"
for mode in fixed scaled; do
    python "$p/c01_score_route.py" "$w/$mode" "$w/c01/$mode" "$mode"
done
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
args=(--edition 2024 -O -Dwarnings -L dependency=target/release/deps
      --extern soil_fem_3d=target/release/libsoil_fem_3d.rlib --extern "serde_json=$serde")
rustc "${args[@]}" "$p2/nc_map.rs" -o "$w/nc-map"
rustc "${args[@]}" "$p/nc_units.rs" -o "$w/nc-units"
rustc "${args[@]}" --cfg roundtrip_units "$p/nc_units.rs" -o "$w/nc-roundtrip"
for ds in 0.02 0.04 0.32; do
    python "$p/nc_branch.py" "$w/inputs" "$w/nc-map" "$w/nc/$ds" "$ds"
done
python "$p/nc_branch.py" "$w/inputs" "$w/nc-units" "$w/nc/unit-0.04" 0.04
python "$p/nc_branch.py" "$w/inputs" "$w/nc-roundtrip" "$w/nc/roundtrip-0.04" 0.04
python "$p/nc_audit.py" "$w/inputs" "$w/nc" > "$w/nc/audit.json"
cargo test --release --locked --test joint_path_error --test joint_law
