#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "$1")
p=interim/2026-11-26-remaining-p6
p5=interim/2026-11-26-remaining-p5
mkdir "$w/probe"
git ls-files -z -- Cargo.toml Cargo.lock rust-toolchain.toml src examples | tar --null -T - -cf - | tar -xf - -C "$w/probe"
cp "$p5/nc_point.rs" "$w/probe/examples/p5_point.rs"
cp interim/2026-11-25-remaining-p4/nc_homogeneous.rs "$w/probe/examples/p5_homogeneous.rs"
(cd "$w/probe" && git apply -p1 - < "$OLDPWD/$p5/observe.patch")
python "$p/instrument.py" "$w/probe"
(cd "$w/probe" && cargo build --locked --release --example c01_methods --example p5_point --example p5_homogeneous)
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
