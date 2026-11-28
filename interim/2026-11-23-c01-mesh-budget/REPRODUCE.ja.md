# C-01 予算比較の再現

リポジトリの固定依存を導入し、ルートで以下を Bash から実行する。新しい出力先を使用する。Windows の個別実行では `--binary target/release/examples/c01_methods.exe` を指定する。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
driver=interim/2026-11-23-c01-method-improvements/run_case.py
study=interim/2026-11-23-c01-mesh-budget
out=work/new-budget-routes
seeds=work/new-budget-seeds
cargo build --locked --release --example c01_methods
for side in native independent; do
  for n in 256 512; do
    python "$driver" "$seeds/$side-$n-base" --solver "$side" --steps "$n" --save-seed
  done
  for mode in base ws; do
    for n in 512 1024; do
      opts=()
      if [[ $mode == ws ]]; then opts=(--active --seed "$seeds/$side-$((n/2))-base"); fi
      python "$driver" "$out/$side-$n-$mode" --solver "$side" --steps "$n" --attempt-budget-policy mesh_proportional --history "${opts[@]}"
    done
    if python -c 'import json,sys; sys.exit(0 if json.load(open(sys.argv[1]))["complete"] else 1)' "$out/$side-1024-$mode/result.json"; then
      opts=()
      if [[ $mode == ws ]]; then
        python "$driver" "$seeds/$side-128-base" --solver "$side" --steps 128 --save-seed
        opts=(--active --seed "$seeds/$side-128-base")
      fi
      python "$driver" "$out/$side-256-$mode" --solver "$side" --steps 256 --attempt-budget-policy mesh_proportional --history "${opts[@]}"
    fi
  done
done
python "$study/measure.py" "$out" "$out/results" --seeds "$seeds"
```

全工程の受理時刻・32 GP の累積塑性散逸を原票に保存する。対照は旧比較の同一実装・N・方式。N1024 完走時は同じ方式の N256 を第3水準として用い、T8 の規約どおり、隣接2水準の差は粗い側の基準座標、次数は3水準共通の N256 座標で比較する。全7工程の実受理値を用い、途中の補間・ゼロ埋めは行わない。停止した経路では完走を前提とした次数を算出しない。

PC 実行はしていない。PC で再現する場合は12 Tet10、直列、BLAS 1スレッド、1経路32 GB・6時間の計画上限。768/6144は対象外。
