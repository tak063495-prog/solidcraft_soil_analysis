# C-01 再現・Windows PC 手順

保存先ラベルは依頼どおり 2026-12-11。PR #13、起点 `9338204990ebb529b811918dc173af7aa9572903` の同じブランチ。マージしない。以下の出力先は毎回未存在の名前を用いる。PC の実計算は今回実行していない。

## 準備（リポジトリルート、PowerShell）

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 9338204990ebb529b811918dc173af7aa9572903 HEAD
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
cargo build --locked --release --example wet_joint_step_probe --example wet_joint_history_probe
```

終了コードが非ゼロならログと途中場を保存し、成功として先へ進めない。Linux/WSL は `.venv/bin/python` と拡張子のない実行ファイルを使う。証拠の SHA は実バイトを検証するため、証拠ファイルは `core.autocrlf=false` で取得する。文書を読む新しい試験は LF/CRLF をテキストとして同じに扱う。新フォルダの `.gitattributes` はJSON・ログ・TSV・ZIP partを `-text` にし、autocrlfでも封印した生バイトを保持する。

## 新オプションの小ケースを新規求解

```powershell
python interim/2026-12-11-c01-global-step-control/study.py --mode init --solver target/release/examples/wet_joint_step_probe.exe --output work/NEW-c01-control
python interim/2026-12-11-c01-global-step-control/study.py --mode native --output work/NEW-c01-control
python interim/2026-12-11-c01-global-step-control/study.py --mode independent --output work/NEW-c01-control
python interim/2026-12-11-c01-global-step-control/study.py --mode restart --steps 256 --output work/NEW-c01-control
python interim/2026-12-11-c01-global-step-control/measure.py work/NEW-c01-control work/NEW-c01-control-measure
```

入力は既存 `examples/joint-water-column-tet10-mesh.json`、`ground-joint-water-mc.json`、`stages-joint-water-cycle.json`。正規化したコピーを出力直下へ保存する。固定54節点・12 Tet10・324材料GP・32接触GP、固定水位0.75 mのCONTINUOUS、MC弾性枝、全7工程を保持する。各工程の基準全体刻み N=16/32/64/128/256 に対し、`joint_system.global_step_control={"switch_bisections":4}` を明示する。接触内上限は5e−5 m、最小全体増分は `min(旧値, 1/(32N))`。これは新入力の経路分解能であり、力1e−8・真の線形1e−10・独立系1e−10を変更しない。

本体・独立参照は入力だけから別々に切替を判断する。全受理場と棄却を保存し、各 `i/N` と工程端点を通る。独立系へ本体解や本体の受理刻み列を渡さない。`gN/native-raw/accepted.jsonl`・`run.json`・`progress.jsonl`・実 `checkpoint-*.json`、可逆NPZの `native/`・`independent/`、監査JSON、再開後の全状態を出力する。最細の7工程からそれぞれ実ファイルを別PIDで復旧し、全後続確定状態とRun全体を時間だけ除いて比較する。`--whole` の中断なし対照も照合する。

想定は一式20〜60分、RAM16 GB、ディスク10〜20 GB。CPU・ディスク・観測用checkpoint処理に依存する資源計画で、Windows実測値ではない。1プロセスの既定予算は1800秒。未収束・時間切れ・試行上限・欠落を成功にしない。見る値は `measurements.json` の全18場の全履歴L2、同一座標の絶対差から求めた観測次数・差方向、statusカテゴリ不一致、`path_audits` の元の残差、`separate_process_restart`。単調になっても自動でC-01を受入にしない。

## 既存記録の --maximum 0.000025 を正確に実行する PC コマンド

次は **旧方式の接触内上限だけを2.5e−5 mにする対照**。新しい全体刻み制御の有効化とは別の計算である。入力は上記と同じ3つの examples であり、ドライバが出力直下の `mesh-h0.json`・`ground.json`・`stages.json` に保存する。

```powershell
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c01-g256-512-1024-m000025 --levels 0 --steps 256 512 1024 --maximum 0.000025 --native-timeout-seconds 14400
python interim/2026-11-18-wet-joint-step-convergence/pc_summary.py --input work/NEW-PC-c01-g256-512-1024-m000025 --output work/NEW-PC-c01-g256-512-1024-m000025-summary --steps 256 512 1024 --h 0 --maximum 0.000025
```

出力先は `work/NEW-PC-c01-g256-512-1024-m000025/h0-g256-m2.5e-05-cut/`、同 `g512`・`g1024` の本体Run・独立参照・全受理場。集計は `work/NEW-PC-c01-g256-512-1024-m000025-summary/measurements.json`。見る値は `native` / `independent` の256→512→1024の塑性すべり・累積すべり・全節点/全材料GPのL2と絶対差、共通座標の観測次数・方向、元の残差、全工程完了・カテゴリ差、`passed=false` の維持。上限5e−5 mの別ディレクトリと比較して局所上限依存を分離する。想定30分〜2時間以上、RAM16 GB、ディスク10〜30 GB、本体1経路の予算4時間。未実測のPC計画であり、保証時間ではない。

## PC のみ：新方式の長い全体系列

```powershell
python interim/2026-12-11-c01-global-step-control/study.py --mode init --solver target/release/examples/wet_joint_step_probe.exe --output work/NEW-PC-c01-control-long
python interim/2026-12-11-c01-global-step-control/study.py --mode native --pc-series --steps 256 512 1024 --timeout 14400 --output work/NEW-PC-c01-control-long
python interim/2026-12-11-c01-global-step-control/study.py --mode independent --pc-series --steps 256 512 1024 --output work/NEW-PC-c01-control-long
python interim/2026-12-11-c01-global-step-control/study.py --mode restart --pc-series --steps 1024 --timeout 14400 --output work/NEW-PC-c01-control-long
python interim/2026-12-11-c01-global-step-control/measure.py work/NEW-PC-c01-control-long work/NEW-PC-c01-control-long-measure --steps 256 512 1024
```

想定1〜4時間以上、RAM16 GB、ディスク30 GB。全体の2000試行/工程という既存予算は保持しているため、細分化により予算が枯渇した場合も未達として保存する。外挿やRichardson補正で合格を作らない。

## PC のみ：Tet10 空間・時空間の生データ

```powershell
python interim/2026-11-15-wet-joint-convergence/user_pc.py --solver target/release/examples/wet_joint_history_probe.exe --output work/NEW-PC-c02-h012-g64-128-256-m000025 --levels 0 1 2 --steps 64 128 256 --maximum 0.000025 --native-timeout-seconds 14400
```

12/96/768 Tet10 と全体64/128/256の交差系列。出力は `mesh-h0.json`〜`mesh-h2.json` と `hH-gN-m2.5e-05-cut/`、最上位 `measurements.json`。RAM64 GB以上、ディスク100 GB、数時間〜1日を計画する。既存ドライバは親面平均と共通節点の診断であり、C-02の共通物理点・全場の空間受入の代用にはならない。C-02〜C-04、NC-01、P6、A/B各項は未達を維持する。

## 検査と台帳

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

Pythonはtests/toolsの全発見IDと4群の実行IDの和集合が一致することを検査する。PowerShell試験は実物の `pwsh` をPATHへ用意する。skipで未導入を回避しない。保存証拠の展開もあるため全件は数十分〜1時間以上とディスク余裕が必要。件数・失敗・未実行・skipは `validation.json` を見る。push後のWindows CIは未確定で、Linux成功や旧先端成功を転用しない。

新しい `acceptance.json`・`measurements.json`・`validation.json` は、先に不変コミットへ保存してから `interim/2026-12-10-ci-p601-inventory-fix/register_sources.py --commit FULL_SHA --id UNIQUE_ID interim/2026-12-11-c01-global-step-control/FILE.json` で生成した登録を同台帳へ追記し、`SUPPLEMENTAL_SHA256` を更新する。旧登録・旧JSONは再生成しない。登録後はP6-01/P6-03監査を実行する。生バイトのSHA検証と文書のLF/CRLF比較を混同しない。

C-01は未達、passed=false。
