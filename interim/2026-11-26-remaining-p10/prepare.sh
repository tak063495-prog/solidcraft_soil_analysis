#!/usr/bin/env bash
set -euo pipefail
w=$(realpath "${1:?fresh scratch}")
p=interim/2026-11-26-remaining-p10
p5=interim/2026-11-26-remaining-p5
mkdir "$w/probe"
git ls-files -z -- Cargo.toml Cargo.lock rust-toolchain.toml src examples | tar --null -T - -cf - | tar -xf - -C "$w/probe"
cp "$p5/nc_point.rs" "$w/probe/examples/p5_point.rs"
cp interim/2026-11-25-remaining-p4/nc_homogeneous.rs "$w/probe/examples/p5_homogeneous.rs"
(cd "$w/probe" && git apply -p1 - < "$OLDPWD/$p5/observe.patch")
python interim/2026-11-26-remaining-p6/instrument.py "$w/probe"
cp src/stages_joint_path_error.rs "$w/probe/src/"
python interim/2026-11-26-remaining-p7/prepare.py "$w/probe"
python interim/2026-11-26-remaining-p8/prepare_contact.py "$w/probe"
python interim/2026-11-26-remaining-p9/instrument.py "$w/probe"
python "$p/instrument.py" "$w/probe"
(cd "$w/probe" && cargo build --release --locked --example c01_methods --example p5_point --example p5_homogeneous)
python interim/2026-11-22-impl-linear-solver/prepare_inputs.py "$w/inputs" --cloud-original
