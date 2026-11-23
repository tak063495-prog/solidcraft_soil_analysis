# 範囲別判定と原因調査の再現

リポジトリのルートで実行する。固定 Python 環境と lock 済み Rust・NumPy・SciPy を使用する。出力先は未使用のディレクトリを指定する。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
D=interim/2026-11-23-ac-scope-decisions
C=interim/2026-11-23-impl-root-causes-c
W=runs/ac-decisions-replay
mkdir -p "$W"
cargo build --locked --release --example t7_boundary --example joint_a02_projection_probe \
  --example c03_bulk_cause_probe --example c03_bulk_return_trace
python "$D/open_trials.py" target/release/examples/t7_boundary "$W/a01.json"
python "$D/shear_transfer.py" target/release/examples/joint_a02_projection_probe "$W/a02"
M="$W/matrix"
python "$C/c03_cloud.py" "$M" --initialize target/release/examples/c03_bulk_cause_probe
for side in native independent; do
  for route in m1c1 m0c1 m1c0; do
    python "$C/c03_cloud.py" "$M" --route "$route" --solver "$side"
  done
done
python "$D/native_sensitivity.py" "$M"
python "$D/classifications.py" "$M" target/release/examples/c03_bulk_return_trace "$W/c03.json"
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root .
```

A-01 は固定 ZIP の SHA-256 を再照合し、両系の保存親状態から各候補を再評価する。予測点の評価には新しい Newton 試行を追加しない。GP 番号は独立の物理順で0始まり、本体は所有要素対・面・物理位置・法線・重みによる一対一対応を使用する。出力の `stored_label_mismatches` は保存候補の実ラベルと追跡器の到着ラベルの相違を示す。

A-02 の request.json は二つの載荷面の接触履歴から計算した節点力。native.json に本体の受理状態・GP配置・受け渡し後の全節点力、comparison.json に解析積分との比較を記録する。

C-03 は 12 Tet10 に限定する。元の本体3経路を診断出力付きで再実行し、全受理状態のバイトが通常実行と一致することを確認する。独立側は各候補を独立に求解し、カットバックがあれば直前の追加受理点を親に使う。本体の厳しい残差条件での3経路は別フォルダに保存する。`classifications.py` は7組の全224共通座標の分類差が旧照合と同じかを検査する。局所交換では既存の本体材料プローブと独立材料関数を用い、物理受入の履歴に交換結果を加えない。
