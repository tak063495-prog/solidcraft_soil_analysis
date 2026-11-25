#!/usr/bin/env bash
# Run from repository root in the locked toolchain. Only the scratch copy is patched.
set -euo pipefail
p=interim/2026-11-26-remaining-p5
p4=interim/2026-11-25-remaining-p4
w=$(realpath "$1")
mkdir "$w/probe"
files=(Cargo.toml Cargo.lock rust-toolchain.toml src examples)
git diff --quiet 52648fce72dadee2f2bfb89c8cca0f108cb4d3b9 -- "${files[@]}"
git ls-files -z -- "${files[@]}" | tar --null -T - -cf - | tar -xf - -C "$w/probe"
cp "$p/nc_point.rs" "$w/probe/examples/p5_point.rs"
cp "$p4/nc_homogeneous.rs" "$w/probe/examples/p5_homogeneous.rs"
cp interim/2026-11-25-remaining-p3/nc_units.rs "$w/probe/examples/p5_units.rs"
(cd "$w/probe" && git apply -p1 - < "$OLDPWD/$p/observe.patch" && cargo build --locked --release --example c01_methods --example p5_point --example p5_homogeneous)
(cd "$w/probe" && cargo rustc --locked --release --example p5_units -- --check-cfg 'cfg(roundtrip_units)' --cfg roundtrip_units)
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
