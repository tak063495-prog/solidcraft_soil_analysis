# C-01 新規求解・再現とWindows PC手順

PR #13、`codex/plastic-spatial-convergence`、起点 `2219f906431b66927fdbb39864372f2e8e2ad22d`。保存先ラベルは依頼どおり2026-12-12。マージしない。各出力先は未存在の名前を用いる。PCの計算は今回実行していない。

## 準備（リポジトリルート、PowerShell）

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 2219f906431b66927fdbb39864372f2e8e2ad22d HEAD
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements-reference.lock.txt
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
cargo build --locked --release
cargo build --locked --release --example wet_joint_path_probe --example wet_joint_step_probe --example wet_joint_history_probe
```

終了コード非ゼロはそのまま失敗として保存する。Linux/WSLでは `.venv/bin/python` と拡張子なしの実行版を使う。証拠のSHAは実バイトを照合する。新フォルダの `.gitattributes` はJSON/JSONL/ログ/TSV/partを `-text` とし、文書の比較は `read_text` の改行正規化でCRLFも扱う。証拠SHAを改行正規化で弱めない。

## 同じ小ケースの全体16/32/64/128/256

```powershell
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver target/release/examples/wet_joint_path_probe.exe --output work/NEW-c01-path-error
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --output work/NEW-c01-path-error
python interim/2026-12-12-c01-switch-event-localization/study.py --mode independent --output work/NEW-c01-path-error
python interim/2026-12-12-c01-switch-event-localization/restart_one.py --steps 256 --output work/NEW-c01-path-error
python interim/2026-12-12-c01-switch-event-localization/measure.py work/NEW-c01-path-error work/NEW-c01-path-error-measure
```

元入力は `examples/joint-water-column-tet10-mesh.json`、`ground-joint-water-mc.json`、`stages-joint-water-cycle.json`。54節点・12 Tet10・324材料GP・32接触GP、固定水位0.75 m、CONTINUOUS、MC弾性枝、全7工程を保持する。新オプションは次の明示入力だけで有効になる。旧 `switch_bisections` との同時指定は拒否する。

```json
{"global_step_control":{"path_error":{"slip_scale_m":0.0001,"dissipation_scale_kpa_m":0.01}}}
```

接触内上限5e-5 m、最小全体増分 `min(旧値,1/(65536N))`。scaleは積分精度の数値パラメータであり、受入許容値ではない。各GPの1増分対2半増分の差を、すべりの成分最大値と散逸の絶対値で測り、`base_increment²*scale` を超えれば全体増分を棄却する。受理後の提案刻みは誤差比から変更するが、毎回同じ条件で再判定する。元の力1e-8・真の線形1e-10・独立系1e-10、2000試行/工程の上限を保持する。局所jumpや保存解の内分、外挿・Richardson補正は使わない。

出力は `gN/native-raw/` の確定全状態JSONL、全Run、実checkpoint、進捗、`gN/native/`・`independent/` の可逆NPZ、監査JSON、独立試行の `independent-attempts.jsonl`。数値停止時は `failed-run.json` または `independent-failure.json` と途中場を残し、完了やL2を捏造しない。見る値は全工程の完了、元の残差、誤差推定の成分/GP・budget・途中カテゴリ・棄却、全18場の全履歴L2、3水準の観測次数・差方向余弦、切替カテゴリ差、実再開の一致である。旧系列・2026-12-11系列は保存値として並記し、書き換えない。

想定20〜90分、RAM16 GB、ディスク10〜30 GB（履歴・checkpoint数に依存）。Windowsでの実測時間ではない。1経路の既定予算1800秒。中断・未収束・最小刻み・試行上限は失敗のまま記録する。

## 保存証拠の再計測

```powershell
python interim/2026-12-12-c01-switch-event-localization/evidence.py --output work/NEW-c01-event-replay
```

part・ZIP・全エントリーのSHAを確認して展開し、保存全場から計測JSONと切替イベントを再計算する。これは新規数値求解やPC実行ではない。実checkpoint再開は上の `restart_one.py` で最細経路の実ファイル1件から別プロセスで完了まで行う。全7工程からの再開を確認したとはしない。Rust回帰では別途、中間保存後と中断なしの全Runも照合する。

## 旧既定・2026-12-11切替制御のバイト回帰

旧実行版を指定起点のソースだけから別にビルドする。固定参照値・旧checkpointの実行版SHAを付け替えない。

```powershell
New-Item -ItemType Directory work/NEW-base-source
git archive --format=tar --output=work/NEW-base-source.tar 2219f906431b66927fdbb39864372f2e8e2ad22d Cargo.toml Cargo.lock rust-toolchain.toml THIRD_PARTY_NOTICES.md src examples
tar -xf work/NEW-base-source.tar -C work/NEW-base-source
cargo build --locked --release --manifest-path work/NEW-base-source/Cargo.toml --example wet_joint_step_probe
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver work/NEW-base-source/target/release/examples/wet_joint_step_probe.exe --output work/NEW-default-base
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --legacy --steps 16 --output work/NEW-default-base
python interim/2026-12-12-c01-switch-event-localization/compatibility.py --solver target/release/examples/wet_joint_step_probe.exe --baseline-default-raw work/NEW-default-base/g16/native-raw --output work/NEW-legacy-compatibility
```

`compatibility.json` と `default/g16/raw/`、`switch/gN/raw/` が出力先。旧既定112状態の実バイト、切替制御全5水準の確定状態SHAを旧封印索引に照合し、旧12 partの全SHAと全Run（実測時間・実行版SHAだけ除く）も照合する。想定10〜30分、RAM16 GB、追加ディスク3〜6 GB。Windows固有の数値差が出ても許容値で吸収せず、その失敗を残す。Linuxでの今回の結果をWindows成功へ転用しない。

## PCのみ：新方式の長い全体系列

```powershell
python interim/2026-12-12-c01-switch-event-localization/study.py --mode init --solver target/release/examples/wet_joint_path_probe.exe --output work/NEW-PC-c01-error-long
python interim/2026-12-12-c01-switch-event-localization/study.py --mode native --pc-series --steps 256 512 1024 --timeout 14400 --output work/NEW-PC-c01-error-long
python interim/2026-12-12-c01-switch-event-localization/study.py --mode independent --pc-series --steps 256 512 1024 --output work/NEW-PC-c01-error-long
python interim/2026-12-12-c01-switch-event-localization/study.py --mode restart --pc-series --steps 1024 --timeout 14400 --output work/NEW-PC-c01-error-long
python interim/2026-12-12-c01-switch-event-localization/measure.py work/NEW-PC-c01-error-long work/NEW-PC-c01-error-long-measure --steps 256 512 1024
```

想定1〜6時間以上、RAM16 GB、ディスク30〜60 GB。出力は `g256/`・`g512/`・`g1024/` と別集計先の `measurements.json`。全18場の単調性・絶対差・次数/方向、元の残差、カテゴリ、全工程完了・再開一致を見る。2000試行/工程は保持するため、予算停止も未達として残す。単調な塑性すべりだけでC-01を受入にしない。

## PCのみ：接触内上限対照・Tet10空間系列

次は旧方式の接触内上限を変える別の対照である。新しい全体誤差制御の有効化と区別する。

```powershell
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c01-g256-512-1024-m000025 --levels 0 --steps 256 512 1024 --maximum 0.000025 --native-timeout-seconds 14400
python interim/2026-11-18-wet-joint-step-convergence/pc_summary.py --input work/NEW-PC-c01-g256-512-1024-m000025 --output work/NEW-PC-c01-g256-512-1024-m000025-summary --steps 256 512 1024 --h 0 --maximum 0.000025
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c02-h012-g64-128-256-m000025 --levels 0 1 2 --steps 64 128 256 --maximum 0.000025 --native-timeout-seconds 14400
```

前二つは固定12 Tet10、`h0-gN-m2.5e-05-cut/` の本体・独立場と別集計先のL2/次数/方向を、上限5e-5の保存系列と比較する。想定30分〜2時間以上、RAM16 GB、ディスク30 GB。最後は12/96/768 Tet10×全体64/128/256。`mesh-hH.json` と `hH-gN-m2.5e-05-cut/` が出力先、想定数時間〜1日、RAM64 GB以上・ディスク100 GB。親面平均や共通節点の診断はC-02の共通物理点・全場受入の代用にはならない。大規模系列はクラウドでもPCでも今回実行していない。

## 全件検査・台帳

```powershell
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-c01-check/core0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-c01-check/core1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-c01-check/mean
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-c01-check/spatial
```

tests/toolsの全発見IDと4群の実行IDの和集合を照合する。実物の `pwsh` が必要で、未導入をskipにしない。全件は数十分〜1時間以上、保存証拠を含むcheckoutに約26 GB、Gitオブジェクトと試験一時領域にも十分な空きを用意する。成功/失敗/未実行/skip・ログ・前後指紋は `validation.json` を見る。push後の新先端Windows CIは未確定。

新しい状態証拠JSONは先に不変コミットへ保存し、`interim/2026-12-10-ci-p601-inventory-fix/register_sources.py --commit FULL_SHA --id UNIQUE_ID interim/2026-12-12-c01-switch-event-localization/FILE.json` の出力を既存追補台帳に追加する。旧登録は再生成せず、`SUPPLEMENTAL_SHA256` を更新してP6-01/P6-03を再検査する。SOURCE_NOT_INVENTORIEDの例外やskipは追加しない。

C-01は未達、passed=false。
