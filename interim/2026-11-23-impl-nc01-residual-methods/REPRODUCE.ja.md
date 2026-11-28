# NC-01 の96 Tet10を同じ実行版で再計算する

Windows の PowerShell で、PR #13 の作業コピーを使う。NC-01 の元入力3ファイルと過去の出力を保存し、6経路とも新しいフォルダで初期状態から開始する。PC 上の計測値はクラウドと別に保存する。

## 実行版と入力

Rust は rust-toolchain.toml の 1.99.1、Python は 3.13。実測時の実装は Git bundle から取得し、集計器は原票の固定コミットから読み出す。台帳の照合には自動改行変換を無効にした作業コピーを使う。以下の最初のコマンドは未使用のディレクトリ名で実行する。

```powershell
git clone -c core.autocrlf=false --filter=blob:none --branch codex/plastic-spatial-convergence https://github.com/tak063495-prog/solidcraft_soil_analysis.git nc01-methods
if ($LASTEXITCODE) { throw 'clone failed' }
Set-Location nc01-methods
$T = 'interim/2026-11-23-impl-nc01-residual-methods'
$evidence = '6ae41048426c86a56708ed246d22160e610d9930'
$implementation = 'c7397678e64a2899a0097fb54ec17d45c0e19bef'
$bundle = "$T/measured-source.bundle"
if ((Get-FileHash -Algorithm SHA256 $bundle).Hash -ne 'd2b1631be12dafe128705dca46b0ec2346820c2b48f24828e39c06667910e2a3') { throw 'bundle hash mismatch' }
git bundle verify $bundle
if ($LASTEXITCODE) { throw 'bundle or base commit missing' }
git fetch $bundle refs/heads/nc01-measured
if ($LASTEXITCODE -or (git rev-parse FETCH_HEAD) -ne $implementation) { throw 'measured commit missing' }
git switch --detach $implementation
if ($LASTEXITCODE -or (git rev-parse HEAD:src) -ne '832a3a6c711b61c268d203a75efc37f6ad49b2d9') { throw 'measured source mismatch' }
git restore --source=$evidence --worktree -- "$T/summarize.py"
if ($LASTEXITCODE) { throw 'summary helper missing' }
git rev-parse HEAD
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE) { throw 'dependencies failed' }
cargo build --locked --release --example nc_history_trace
if ($LASTEXITCODE) { throw 'build failed' }
$Trace = (Resolve-Path 'target/release/examples/nc_history_trace.exe').Path
Get-FileHash -Algorithm SHA256 $Trace
$Out = 'work/nc96-methods-01'
python "$T/prepare.py" $Out --cloud-original
if ($LASTEXITCODE) { throw 'input preparation failed' }
& "$T/run-pc.ps1" -InputDirectory $Out -Trace $Trace
```

PC に保存した元入力を使う場合は、新規の出力名で `--cloud-original` を `--source 'D:/pc_runs/fine_20261122/nc96/input'` に置き換える。このフォルダには初期の mesh.json、ground.json、stages.json が必要。生成器は3ファイルをバイトのまま保存して SHA-256 を表示する。元入力の3 SHA、Git SHA、実行ファイルの SHA を実行ごとに残す。

各経路の `run.json` の `last_committed.stage_progress` が到達 λ、`attempts[-1].progress_to` が最後の目標、`failure_reason` が停止理由。`linear_solves` の `true_relative_residual` を確認する。終了コード0でも `complete=false` の場合がある。棄却試行内の途中の線形報告が Run に残らないことがある。合否には元の binary64 行列積を使い、スケール後の値や高精度観測値を代入しない。

## 予算

96 Tet10・250節点の1経路ずつを処理する。入力の自由度収容数750で密行列1枚は4.5 MB、撤去終端348自由度では0.97 MB。行列と LU を持つ場合の2枚合計は最大9 MB。既往のクラウド96要素本計算は約9～12分。PC の予定は1経路1時間・1 GBを目安とし、スクリプトは6時間または32 GB（ワーキングセットまたはプライベートメモリ）で終了して部分結果を保存する。監視間隔は250 ms。指定上限は1経路1日・40 GB。対象メッシュは96のみ。

## 保存した8系と検証

```powershell
cargo test --locked --release --lib stages::dense_residual -- --nocapture
if ($LASTEXITCODE) { throw 'replay failed' }
cargo fmt --all --check
if ($LASTEXITCODE) { throw 'fmt failed' }
cargo clippy --locked --all-targets -- -D warnings
if ($LASTEXITCODE) { throw 'clippy failed' }
cargo test --locked --release --all-targets --no-fail-fast
if ($LASTEXITCODE) { throw 'Rust tests failed' }
cargo test --locked --release --doc
if ($LASTEXITCODE) { throw 'Rust doc tests failed' }
```

PC の6経路と保存系を同じ形式の新しい JSON に集計するには、同じディレクトリで次を実行する。PowerShell の出力リダイレクトによる UTF-16 化を避け、Cargo の出力をそのまま保存する。

```powershell
python -c "import subprocess,sys; f=open(sys.argv[1], 'xb'); subprocess.run(['cargo','test','--locked','--release','--lib','stages::dense_residual','--','--nocapture'], stdout=f, stderr=subprocess.STDOUT, check=True)" "$Out/replay.log"
if ($LASTEXITCODE) { throw 'saved system replay failed' }
python "$T/summarize.py" $Out "$Out/replay.log" "$Out/pc-measurements.json"
if ($LASTEXITCODE) { throw 'comparison summary failed' }
```

保存系の `NC01_REPLAY` 行には元ファイル、行番号（0始まり）、設定、判定残差、実時間が出る。反復回数は通常 LU 求解を1として数える。P6-01 登録後の PR 先端を新たにチェックアウトした検証用コピーでは、次で台帳と Python 全群を検査できる。

```powershell
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE) { throw 'P6-01 failed' }
python interim/2026-12-07-p6-03-cloud-reconciliation/audit_status.py
if ($LASTEXITCODE) { throw 'P6-03 failed' }
cargo build --locked --release --bins --example nc_input_roundtrip --example nc_dense_lu_diagnostic --example wet_joint_local_path_probe
if ($LASTEXITCODE) { throw 'test binary build failed' }
foreach ($g in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')) {
    python tools/run_ci_python_evidence_shards.py --group $g --output "work/nc96-tests-01/$g"
    if ($LASTEXITCODE) { throw "$g failed" }
}
```

Linux は固定依存に加えて `python -m pip install pytest==9.2.0` を実行する。Linux 専用の7件は `python -m pytest -q tests/linux` で別に実行する。旧文書・台帳の SHA 検証には対応する固定コミットの取得が必要。
