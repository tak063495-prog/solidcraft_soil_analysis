# ユーザーPCでの再計測・新規求解・残件

PR #13の同じブランチを使い、マージしない。作業中の変更を消さず、出力先は毎回未存在のパスにする。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor c71777d2dfd5cde12774981bdd357cf46c69fa06 HEAD
```

最後の終了0は祖先確認だけで、数値受入ではない。今回の環境は [environment.json](environment.json) を参照する。依存は既存のロックファイルを使用する。

```sh
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m pip install -r requirements-reference.lock.txt
```

Linux/WSLでは `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 PYTHONUTF8=1 PYTHONIOENCODING=utf-8` をexportする。PowerShellでは各変数を `$env:OPENBLAS_NUM_THREADS='1'` の形式で設定する。

## 1. 保存証拠の再監査（求解なし）

```sh
python interim/2026-11-10-wet-joint-small-refinement/replay.py --output work/NEW-small-joint-replay
```

期待する記録は121エントリー・789封印ソースのSHA確認、`all_measurements_exact=true`、`accepted_states=1680`、`numerical_acceptance=false`。全1,680保存場を再組立し、隔離UTF-8の別プロセスで封印ソースから計測JSONを再生成する。PC/BLAS差で末尾の算術差が出た場合は、その不一致を記録する。基準・封印参照値を変更して通さない。Windowsでの新規数値再現はこのLinux記録から認定しない。

## 2. 同じ入力から16/32/64/128を新規計算

```sh
python interim/2026-11-10-wet-joint-small-refinement/study.py --output work/NEW-small-joint-solve
python interim/2026-11-10-wet-joint-small-refinement/study.py --replay work/NEW-small-joint-solve --output work/NEW-small-joint-remeasurement.json
```

これは本体場を使わない既存のMC弾性バルク＋Coulomb JOINT参照で、全体分割だけを変える。生産Rustの体積積分や本体checkpointを変更しない。`g16-l8`〜`g128-l8` の `result.json` と `accepted/`、`measurements.json` を残す。16/32の全比較が旧JSONと一致するか、32/64・64/128のすべり相対L2と絶対差が残るかを確認する。今回の値は約0.67627／0.03662／0.07426だが、新しいPCでこの値を入力して結果を作らない。

共通受理点の比較だけでは細経路の全点を比較できない。既存の7端点比較、18場の最大絶対差・相対L2、全点の力・接触伝播・仕事・水・適合ひずみ監査を併記する。`passed=false` と一般塑性・本体比較・完了再開の未確認を維持する。未収束・プロセス異常・時間予算停止の場合は欠損をそのまま残し、7段階完了や物理的破壊にしない。

さらに既存系列を延長する場合は、既存 `bulk.history.cut_volume_reference(*values, steps, 8, 7, history_sink=sink)` の `steps` を256等にし、同じ入力から全前史を計算する。保存は `HistoryWriter`、検証は `validate_path` と `bulk.audit`、128/256比較は `compare_accepted(..., ratio=2)` をそのまま用いる。過去の参照場を初期値や正解として注入しない。接触内分割の影響を調べる場合も既存の8/16等を別対照として分離し、今回の全体分割比較と混ぜない。これらの追加求解は今回実行していない。

## 3. NCと統合の未確認項目

NC全撤去経路の追加求解は今回行っていない。全件試験の既存小ケース求解は別である。既存手順は [材料原因診断の再現手順](../2026-11-10-nc-history-material-cause/REPRODUCE.ja.md)、[既存擬似弧長の再現手順](../2026-11-09-nc-pseudo-arclength/REPRODUCE.ja.md) を参照する。5 kPa追加プレロードを維持する。100 kPaプレロード合格や補助拘束圧を転用しない。

将来の全解放受入には、撤去中の材料・接触履歴、非自明な残存域λ=0のNewton釣合い、元保持荷重ゼロ、全解放後の実ゼロ増分による全履歴凍結、入力のみの独立塑性3D全場、両実装の完了までの刻み感度、実checkpointから完了までの別プロセス再開が必要である。保存ZIPを読み直す検査を、solverの途中再開とは扱わない。checkpointの実行版identityを付け替えたり、保存場を編集・再封印したりしない。

湿潤JOINTも、本体全節点3成分、全材料GP、全接触履歴、物理水荷重と反力・釣合い、同じ荷重経路の独立3D比較が別に必要である。今回のバルクは弾性枝であり、一般バルク塑性や部材との統合を認定できない。証拠不足の合格を拒否する試験は次で確認できる。

```sh
python -m unittest discover -s tests -p test_wet_joint_small_refinement.py -v
python -m unittest discover -s tests -p test_wet_joint_bulk_history.py -v
python -m unittest discover -s tests -p test_wet_joint_accepted_replay.py -v
```

## 4. Linux/Windowsの全件試験

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-check/python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-check/python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-check/python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-check/python-spatial-evidence
```

全件試験には既存mean-stress約5.49 GBとspatial約0.92 GBの分割証拠・展開領域が必要である。今回の33.76 MBのZIPだけでは全件試験の依存を満たさない。non-cone sparse checkoutで省略している場合だけ次を実行する。

```sh
git sparse-checkout add '/interim/2026-11-10-nc-mean-stress-release/*.part*' '/interim/2026-11-10-nc-spatial-release/*.part*'
```

終了0、fail/error/ignored/skipped/filtered-out 0と、各inventoryの発見全ID＝4群の実行ID和集合（重複なし）を確認する。別Windows試験があるためLinuxと件数が異なる。大きい保存証拠のSHA確認・再計測は新規大規模求解ではない。CIの成功は数値全履歴受入と区別する。

## 5. クラウドで実行しない本計算

Griffiths 1200/4800/19200 Tet10、乾燥MCの4/8/16と8/16/32、MCC 274〜14909の再求解、斜面1200要素級SRM、仮想仕事上端の再計算、82944降雨、105456性能、2807要素級実CAD、長時間Richards、実断面SRMの本計算、Native UI、別PC/GPU、完成配布はユーザーPC向けとする。今回これらは実行していない。

既存の各入力・参照値・受入手順を使用し、入力SHA、solver identity、収束・中断理由、全場・履歴、粗密差、実再開、資源計測を保存する。形状や参照値、部分亀裂先端・接触BIRTH・高次板・弱結合の新仕様を推測して補わない。P3b/P4各行、構造・水・JOINT全体、P6全体は、個別に必要な独立3D証拠と試験が揃うまで未完とする。
