# T8 比較の再現

## 起点と環境

PR #13 の履歴、`120a66520b628e1a756e719c374f28d16d9e7581`、登録の固定コミット、各原本アーカイブが必要。以下はリポジトリルートで実行する Linux/WSL の bash。PowerShell 依存テストは 7.5.4 を使用する。Windows native では solver に .exe を付け、UTF-8 のテキスト読取りで改行を正規化する。SHA 検査は生バイトで行う。過去診断の期待 SHA はそれを保存したコミットで検証する。

```sh
set -eu
git fetch origin refs/heads/codex/plastic-spatial-convergence
git merge-base --is-ancestor 120a66520b628e1a756e719c374f28d16d9e7581 HEAD
rustc --version
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export RUST_TEST_THREADS=2 CARGO_BUILD_JOBS=2
T=interim/2026-11-22-impl-t8-acceptance-criteria
W="$PWD/work/NEW-T8"
mkdir -p "$W"
cli="$PWD/target/release/soil-fem-3d"
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bins --examples
```

科学計算には未使用の出力先を使う。確定済みの出力には追記・上書きしない。C-01 は最初に N=128 の経過時間とピーク RSS を見て次水準に進む。今回の Linux 環境のメモリ制限は16 GiB。許容規模は RAM 40 GB 以下・1経路1日以内。768/6144 のソルバー実行はない。旧アーカイブの読取りには過去の生データを含む。C-02 の抽出器は全 part の SHA を検証する。新計測 `two_mesh.py` は h0/h1 のみを使う。

## 封印された比較だけを再検査

```sh
python "$T/observations.py"
python -m unittest discover -s tests -p test_t8_acceptance_criteria.py -v
```

封印には比較の入力になった小型数値結果・拒否記録・刻みごとの有限参照を収録した。出力バイナリや全経路の新規生チェックポイントは含めていない。下記で再生成できる。A-01 は T7 の SHA 検証済み分割 ZIP の元ケース・両解の受理場、弱空間と区間仕事の記録を読む。NC-01 は起点に含まれる T9 の状態 JSON を読む。どちらも追加の NC 求解を起動しない。

## A・B 群

```sh
python "$T/water.py" "$W/t8-water" "$cli"
python interim/2026-11-21-impl-t3-tip-crossing-birth/study.py "$W/a02"
python interim/2026-11-20-impl-t2-soil-birth-hydraulic-init/study.py "$W/a03"
python interim/2026-11-20-impl-t1-wet-joint-richards/study.py controls "$W/a04"
S=interim/2026-11-21-impl-t5-curved-tri6-shell
mkdir -p "$W/b01"
python "$S/modes.py" target/release/examples/curved_shell_operators "$W/b01/modal"
python "$S/manufactured.py" "$cli" "$W/b01/mms"
python "$S/standard.py" "$cli" "$W/b01/standard"
python "$S/attached.py" "$cli" "$W/b01/attached"
python interim/2026-11-21-impl-t4-offset-plate-composite/study.py "$W/b02-patch"
python "$T/plate_area.py" "$W/b02-area"
J=interim/2026-11-21-impl-t6-offset-kink-joint
python "$J/t6_study.py" "$cli" "$W/b03"
python "$J/t6_convergence.py" "$W/b03"
```

A-03 は水面上下の実体積・保水量・重量の原票を `placement.json` に保存する。A-04 の解析参照は H=2 または H=2−0.2z、q=−K∇H、面上の一定圧 γ(H−z) と三角形の各節点荷重 pA/3。両面で外向き流量の符号が反対になる。両試料とも流体供給源0で、機械変形から水理への非ゼロ逆伝達は含まない。各物理時刻の符号を保った差と従来の両水収支を `affine.json` に保存する。

B-02 は同じ荷重・支持・腕の `beam-positive` と `mesh-1` の保存済み本体・独立場を展開する。板各三角形の6点求積で N、M、M−eN、Q を全体軸で面積積分する。接続点からの物理距離は `probe_distance_m`、節点反力と外力仕事の各時点の値は `comparison.json` にある。

## C 群

```sh
S=interim/2026-12-11-c01-global-step-control/study.py
python "$S" --mode init --output "$W/c01" --solver target/release/examples/wet_joint_step_probe
for n in 128 256 512; do
  python "$T/time_case.py" "$W/c01/native-$n.time" python "$S" --mode native --steps "$n" --pc-series --timeout 86400 --output "$W/c01"
done
python "$T/load_steps.py" "$W/c01" "$W/c01-results.json"
python interim/2026-12-17-c02-spatial-cause-diagnosis/c02diag_extract.py --output "$W/c02-source"
python "$T/two_mesh.py" "$W/c02-source" "$W/c02-extended.json"
D=interim/2026-12-18-c03-plastic-bulk-cause-diagnosis
mkdir -p "$W/c03"
python "$T/time_case.py" "$W/c03/matrix.time" python "$D/c03diag_driver.py" --mode matrix --output "$W/c03/matrix" --solver target/release/examples/c03_bulk_cause_probe --element tet4 --steps 32 --timeout 7200
python "$D/c03diag_scope.py" --matrix "$W/c03/matrix" --output "$W/c03/scope" --cli "$cli"
for label in mcc-pc200 mcc-pc105; do
  python "$D/c03diag_mcc_driver.py" --case "$W/c03/scope/$label" --mesh "$W/c03/matrix/mesh.json"
done
python "$D/c03diag_mcc_open.py" --matrix "$W/c03/matrix" --source "$W/c03/scope/mcc-pc105" --output "$W/c03/scope/mcc-pc105-open1" --cli "$cli"
python "$D/c03diag_mcc_driver.py" --case "$W/c03/scope/mcc-pc105-open1" --mesh "$W/c03/matrix/mesh.json"
python "$D/c03diag_measure.py" --matrix "$W/c03/matrix" --scope "$W/c03/scope" --output "$W/c03/results.json"
python "$D/c03diag_replay.py" --matrix "$W/c03/matrix" --output "$W/c03/local-replay"
D=interim/2026-12-19-c04-transient-water-cause-diagnosis/c04cause_driver.py
python "$D" --initialize --output "$W/c04" --probe target/release/examples/wet_joint_transient_probe --solver "$cli"
for pair in '32 8' '32 16' '32 32' '16 32' '64 32'; do
  read -r m w <<< "$pair"
  for side in native adaptive; do
    label="$side-m$m-w$w"
    python "$T/time_case.py" "$W/c04/$label.time" python "$D" --output "$W/c04" --label "$label" --mechanical-steps "$m" --water-steps "$w" --timeout 7200
  done
done
python "$T/contact_labels.py" "$W/c04" "$W/c04-results.json"
mkdir "$W/report"
python "$T/assess.py" "$W" "$W/report"
python "$T/observations.py" --source "$W/report" --pack "$W"
```

C-01 の独立128/256は旧 switch 封印、512は独立原因診断の switch-n512-q4。各 chunk の SHA を検査してデコードし、元の7工程・全32 GP の実在する i/N のみを比較する。試行上限の旧 path-error 512 は別経路。C-03 の新規独立全体求解は本体変位・局所分割数を入力に使わない。凍結経路の再生は `local-replay` に隔離した。C-04 の `boundary` は各不一致の工程・進捗・GP・終点ラベル・両実装の摩擦円余裕・直前小刻みの塑性増分を保持する。原本の不一致数は `preserved_mismatch_counts`、新規の厳密混同行列は `confusion`。旧3%案は両者の判定に使わない。

## 状態監査と全試験

追補台帳は [e541fe1 の登録方式](../2026-12-10-ci-p601-inventory-fix/README.ja.md) に従う。先行コミットに固定した acceptance.json の SHA-256・Git blob・型付き状態を照合し、既存登録の末尾に追加する。監査器には台帳全体の SHA-256 を固定する。再現時には既存の登録を追加し直す必要はない。

```sh
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root .
cargo test --locked --release --all-targets --no-fail-fast
cargo test --locked --release --doc
for group in python-core-0 python-core-1 python-c01-switch python-c01-path python-c01-cause python-spatial-evidence python-mean-evidence; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "$W/$group"
done
python -m pytest -q tests/linux
```

保存した実行記録は `verification.zip` の `summary.json` に従って SHA-256 を検査する。各群の台帳で、その群に割り当てた試験IDと実行時間のある試験IDが一致すること、7群の和集合が core の最終発見IDを一度ずつ含むことを照合する。

Linux の mean/spatial で許容されている既知の setUpClass エラーが出たときはログ・未実行件数を記録して次の群を実行する。成功件数に加えない。これは数値条件の例外ではない。

A-01～C-04、NC-01、P6 の元項目は未達。
