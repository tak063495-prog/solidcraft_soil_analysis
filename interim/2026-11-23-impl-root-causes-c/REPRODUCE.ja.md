# C-01・C-02・C-03 の比較手順

PR #13 の対象ブランチを使用する。リポジトリのルートで以下を実行する。すべて明示選択の診断で、出力先には新しいディレクトリを用いる。新規求解するメッシュは12 Tet10 のみ。C-02 は保存結果を読む。

```sh
T=interim/2026-11-23-impl-root-causes-c
W=work/new-c-comparisons
mkdir -p "$W/trials" "$W/results"
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
cargo build --locked --release --bins --examples
```

## C-01　ギャップ交差と局所積分

各経路を初期平衡から解く。N は時間方向の基準分割数で、要素数ではない。工程上限は各2,000試行。基準は event4、棄却後だけ二分する対照は bisect4、最少内分割だけ1にする対照は event1。どの経路も受理後の基準点への着地規則は同じ。本体は reload の終端で観測を打ち切る。独立系は元手順で全7工程または工程上限まで進む。比較に使う工程は open、reload。出力の `first_rejected_after_final` が null の場合、最終受理後に実際に棄却された目標は記録されていない。

```sh
for n in 512 1024; do
  for route in event4 bisect4 event1; do
    python "$T/c01.py" --steps "$n" --route "$route" "$W/trials/$n-$route.json.gz" || exit 1
    python "$T/c01_native.py" --steps "$n" --route "$route" \
      --solver target/release/examples/joint_path_window --work "$W/native-$n-$route" \
      "$W/trials/native-$n-$route.json.gz" || exit 1
  done
done
python "$T/c01_report.py" "$W/trials" "$W/results/c01.json"
```

同じ受理座標までに費やした棄却には、その座標を一度行き過ぎた試行も含む。試行原票の `gap_zero` は旧状態・中点・細分終点の符号が変わった半区間内の最初の線形予測座標と GP（0始まり）。記録した予測座標は受理状態への補間に使わない。`crossing` はギャップ交差または旧・中点・粗終点・細終点の接触分類の相違、`smooth` はこれらの標本で切替なし、`unknown` は局所誤差観測のない失敗を表す。すべりは各 GP の成分絶対値の最大（m）、散逸は各 GP の絶対差の最大（kPa·m）。受理長は工程進捗の差。

## C-02　同じ448時点・160点の読取り

```sh
D=interim/2026-12-17-c02-spatial-cause-diagnosis
python "$D/c02diag_extract.py" --output "$W/c02-source"
python "$T/c02.py" "$W/c02-source" "$W/results/c02.json"
```

抽出器は旧86分割 ZIP 全体と選択した各メンバーの SHA-256 を確認する。比較器が読むのは h0/h1 の12/96 Tet10・N64/128のみ。旧アーカイブに含まれる他水準の求解は行わない。各工程について実在する受理時刻 i/64 を厳密に探す。時刻の重複・欠落・逆順はエラー。二格子の求積点の和集合、両格子の面積重みを半分ずつ、同一の親面 ID と同一の三次再構成を4対に使用する。点差、面積重み付きの点差、各時点で符号を保って面積積分してから取る差は別欄。最大点差の親面・位置・時刻、および最大積分差時点の左右の積分と符号付き成分差を JSON に保存する。

## C-03　材料／接触細分と不一致 GP

```sh
M="$W/c03-matrix"
python "$T/c03_cloud.py" "$M" --initialize target/release/examples/c03_bulk_cause_probe
for side in native independent; do
  for route in m0c1 m1c0 m1c1; do
    python "$T/c03_cloud.py" "$M" --route "$route" --solver "$side" || exit 1
  done
done
map="$M/m1c1-native/raw/run.json"
for side in native independent; do
  python "$T/c03.py" "$M" "m0c1-$side" "m1c1-$side" --axis material \
    --native-map "$map" --output "$W/results/c03-material-$side.json" || exit 1
  python "$T/c03.py" "$M" "m1c0-$side" "m1c1-$side" --axis contact \
    --native-map "$map" --output "$W/results/c03-contact-$side.json" || exit 1
done
for route in m0c1 m1c0 m1c1; do
  python "$T/c03.py" "$M" "$route-native" "$route-independent" --axis solver \
    --native-map "$map" --output "$W/results/c03-$route-solvers.json" || exit 1
done
```

材料上限は m0=1e−4、m1=5e−5、接触上限は c0=1e−4 m、c1=5e−5 m。全体刻みは32、7工程、open 係数0.12。元の m0c0 カットバック完走を再実行せず、m1c1 を共有する二対を解く。独立参照の入力は mesh・ground・stages のみで、途中の受理点は substeps.jsonl に保存する。左右とも元の224座標を受理したことを確かめてから比較する。節点変位を当該 GP の形状関数で評価する。接触ラベルには両側の変位と同じ GP の traction、両側要素の材料 GP の stress を対応させる。材料ラベルには同じ GP の stress・変位、その要素を所有者とする接触 GP の traction を付ける。隣接接触面がなければ空配列である。

## PC で既存結果を照合する場合

この作業で PC 実行は行っていない。原票が揃う場合は新たに元の m0c0 を求解せず、同じ mesh・ground・stages の SHA を照合し、同一フォルダ構成で C-03 の `--axis solver` に既存の native・independent を渡す。Python は64ビット、BLASは1スレッド。各経路を直列に実施する。C-01 は12 Tet10・32接触 GP でメモリ計画16〜32 GB、時間は1経路6時間。C-03 は同じ12 Tet10・324材料 GP、32 GB、1経路4時間。停止・失敗時の途中記録は保持する。Windows では native 実行ファイルに .exe を付ける。前記の時間は計画上限であり、RAMが40 GBまたは1経路の上限時間に達する前に停止する。C-02 の768/6144は実行しない。[元の資源見積もりと PC 手順](../2026-11-22-impl-c01-c03-nonlinear/REPRODUCE.ja.md)

## 試験と追補台帳

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets --no-fail-fast
cargo test --locked --release --doc
python -m unittest discover -s tests -p test_root_causes_c.py -v
for g in python-core-0 python-core-1 python-c01-switch python-c01-path python-c01-cause python-spatial-evidence python-mean-evidence; do
  python tools/run_ci_python_evidence_shards.py --group "$g" --output "$W/$g"
done
python -m pytest -q tests/linux
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root .
```

Windows の合否は対象コミットの全ジョブ終了時の結果を確認する。台帳は固定コミットから読んだ validation.json の Git blob・SHA-256・型付き状態と、追補台帳全体の SHA-256 を照合する。
