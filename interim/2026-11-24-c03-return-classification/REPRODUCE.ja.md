# C-03 の分類修正前後を再現する

リポジトリのルートで固定 Rust・Python 依存を使用する。出力先は新しいディレクトリを指定する。基準入力は 12 Tet10、m0/m1 は 1e−4/5e−5、c0/c1 は 1e−4/5e−5 m。構造計算の反復上限、荷重刻みと停止条件は基準入力に含まれる。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
D=interim/2026-11-24-c03-return-classification
C=interim/2026-11-23-impl-root-causes-c
S=interim/2026-11-23-ac-scope-decisions
F=interim/2026-11-23-c03-mcc-plastic-work
A=interim/2026-11-23-impl-root-causes-a-nc
P=target/release/examples
W=runs/c03-return-comparison
B="$W/before"
N="$W/after"
M="$W/mcc"
R="$W/regressions"
mkdir -p "$N" "$M" "$R" "$W/work" "$W/results"
cargo build --locked --release --bin soil-fem-3d --example c03_bulk_cause_probe \
  --example c03_bulk_return_trace --example c03_mcc_plastic_work \
  --example t7_boundary --example joint_a02_projection_probe
python "$C/c03_cloud.py" "$B" --initialize "$P/c03_bulk_cause_probe"
cp "$B"/{mesh,ground,stages,identity}.json "$N/"
python "$F/cases.py" "$M" prepare "$P/c03_bulk_cause_probe"
for c in m0c1 m1c0 m1c1; do
  python "$C/c03_cloud.py" "$B" --route "$c" --solver native
  "$P/c03_bulk_cause_probe" "$B/mesh.json" "$B/ground.json" "$B/$c-native/stages.json" "$B/$c-trace" --diagnose
  cmp "$B/$c-native/raw/accepted.jsonl" "$B/$c-trace/accepted.jsonl"
  python "$C/c03_cloud.py" "$B" --route "$c" --solver independent
  python "$C/c03_cloud.py" "$N" --route "$c" --solver independent --resolve-edges
  python "$D/work.py" "$B" "$c" "$W/work/$c-native.jsonl" --native-probe "$P/c03_bulk_return_trace"
  python "$D/work.py" "$B" "$c" "$W/work/$c-before.jsonl" --independent original
  python "$D/work.py" "$N" "$c" "$W/work/$c-after.jsonl" --independent resolved
  python "$F/cases.py" "$M" solve "$c" native
  "$P/c03_mcc_plastic_work" "$M/basis/mesh.json" "$M/$c/raw" "$M/$c/native-work.jsonl"
  python "$F/cases.py" "$M" solve "$c" independent
done
python "$F/compare_work.py" "$M" "$M/results"
python "$D/local.py" "$B" "$P/c03_bulk_return_trace" "$W/results/local.json"
python "$S/open_trials.py" "$P/t7_boundary" "$R/a01.json"
python "$S/shear_transfer.py" "$P/joint_a02_projection_probe" "$R/shear"
python "$A/a02.py" target/release/soil-fem-3d "$R/compression"
python "$A/a04.py" target/release/soil-fem-3d "$R/saturated-longdouble" --reference longdouble
python interim/2026-12-17-c02-spatial-cause-diagnosis/c02diag_extract.py --output "$R/c02-source"
python "$C/c02.py" "$R/c02-source" "$R/c02.json"
python "$D/compare.py" "$B" "$N" "$W/work" "$M" "$R" "$W/results"
python -m unittest discover -s tests -p test_c03_return_classification.py -v
```

本体 MC の `save_substeps` は観測用プローブの要求フィールド。各内増分の終端応力、塑性ひずみ増分、分類を返す。各受理時刻で本体・独立それぞれの材料終端状態を実際の求解結果と厳密照合してから仕事を積算する。本体に関する前後2列は共通の本体計算を用いる。MCC は元の専用プログラムで6経路を独立に再実行する。

`local.json` の材料則交換・親応力交換は保存した修正前の親から始まる局所試験である。`comparison.json` と `mc-work.csv` の修正後は、修正を有効にした独立系で最初から平衡を解いた値である。
