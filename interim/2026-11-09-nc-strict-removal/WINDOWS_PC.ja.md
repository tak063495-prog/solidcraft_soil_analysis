# Windows PC：NC停止の保存証拠監査と、既存小ケースの再実行

物理判定は `passed=false`、`history_converged=false`。Linuxで266増分・0.9514000415802002を新規再現したが、全撤去1.0ではない。この追加コミットのWindows実行・CIは未確認。過去CI 42489252485の成功、20秒予算で進捗0.75のプローブ、粗刻みのNewton上限による停止を、全撤去や今回の到達証拠に転用しない。

PowerShellでリポジトリルートから実行する。既存の出力・保存点を上書きしない。Pythonは既存の固定依存 `requirements-reference.lock.txt`、Rustは `rust-toolchain.toml` の1.99.1を使う。以下のパスは新しい名前に変更できる。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw 'branch update failed' }
git merge-base --is-ancestor 6429e8366ad7b02ab546cc9a9f7ff321b693962b HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong branch history' }
git show -s --format='%H %P' 6429e8366ad7b02ab546cc9a9f7ff321b693962b
# 親: d2b3dd4da9a1e92bbb6bc27b8d03acfc81f0903d
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$Study = 'interim/2026-11-09-nc-strict-removal'
```

## 推奨：求解せず全生証拠を検査

```powershell
$Audit = 'work/NEW-nc-strict-audit'
if (Test-Path $Audit) { throw 'choose a fresh output directory' }
.venv/Scripts/python.exe "$Study/remeasure.py" --output $Audit
if ($LASTEXITCODE -ne 0) { throw 'archive audit differed or failed; retain actual output' }
Get-Content "$Audit/remeasurement.json"
.venv/Scripts/python.exe "$Study/test_audit.py" --evidence "$Audit/evidence" -v
if ($LASTEXITCODE -ne 0) { throw 'audit guard tests failed' }
```

アーカイブ本体と全エントリーのSHA256を検査してから、新規ディレクトリへ展開する。実行バイナリは同梱していない。19＋60＋20線形系をbinary64から正確に取り込み、80桁の独立消去と残差再計測を行う。保存応力の力・保持荷重・履歴・位置保持操作、両実装の部分全場比較も再計測する。Rust・材料写像・独立3D求解を実行する手順ではない。

7つのJSONがLinuxでは完全一致した。WindowsのNumPy/BLAS等の差で一致しなければ、実際の値とエラーを保存する。Linux値を埋めたり閾値を緩めたりしない。`all_matches=true` は再計測の一致であり、`passed=false` と `history_converged=false` は変わらない。

## 任意：未変更本体を新規ビルドし、NCを実際に起動

```powershell
cargo build --locked --release --bin soil-fem-3d
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
$Raw = 'work/NEW-nc-strict-run'
if (Test-Path $Raw) { throw 'choose a fresh raw directory' }
New-Item -ItemType Directory -Path $Raw | Out-Null
.venv/Scripts/python.exe tools/probe_nc_midpoint.py `
  --solver target/release/soil-fem-3d.exe --output "$Raw/fresh-head-midpoint" --budget-seconds 180
if ($LASTEXITCODE -ne 0) { throw 'driver failed; retain all process logs' }
Get-Content "$Raw/fresh-head-midpoint/process-observation.json"
Get-Content "$Raw/fresh-head-midpoint/experiment.json"
.venv/Scripts/python.exe tools/probe_nc_midpoint.py `
  --solver target/release/soil-fem-3d.exe --output "$Raw/fresh-head-resume" --budget-seconds 180 `
  --resume-from "$Raw/fresh-head-midpoint/job"
if ($LASTEXITCODE -ne 0) { throw 'resume observation failed; retain actual evidence' }
Get-Content "$Raw/fresh-head-resume/experiment.json"
```

既存の5 kPa生成器と中点プロファイルを使う。数値差分は `input-controls.json`。物理入力・力1e-8・真の線形1e-10・Newton60を変更しない。Windowsの実測時間・実行版SHA・停止理由・到達進捗を保存する。180秒で途中停止なら、時間切れとして記録する。予算内の終了やドライバ終了コード0を数値収束と扱わない。Linuxでの134.73秒はWindowsの所要時間保証ではない。

通常再開は**そのWindows実行版が作成したジョブ**と同じ封印実行版を使用する。Linuxの封印保存点をWindows実行版へ付け替えたり、元の履歴を修正したりしない。別の未完ジョブから続ける場合も、新しい出力名で `--resume-from` を使用し、実際の再開元を記録する。

## 任意：隔離した最小試行と部分全場比較を再現

```powershell
.venv/Scripts/python.exe "$Study/capture_experiment.py" --output "$Raw/experiment"
if ($LASTEXITCODE -ne 0) { throw 'bounded experiment driver failed; retain actual traces' }
Get-Content "$Raw/experiment/controls.json"
Get-Content "$Raw/experiment/bounded-result.json"
.venv/Scripts/python.exe "$Study/capture_partial_paths.py" `
  --solver target/release/soil-fem-3d.exe --output "$Raw/partial-paths"
if ($LASTEXITCODE -ne 0) { throw 'partial path observation failed' }
.venv/Scripts/python.exe "$Study/audit_evidence.py" --evidence $Raw --output "$Raw/audit-final"
if ($LASTEXITCODE -ne 0) { throw 'measurements differed or incomplete; retain actual failure' }
```

隔離ソースは必ず6429e836から生成される。本体に数値パッチを適用しない。最初に旧ZIP全89エントリーSHAを確認し、保存19系の独立線形ゲートを満たしてから限定試行を行う。古いLinux保存点に対して許される操作は、型付きハッシュと旧識別子を検査する**読み取り専用の限定診断**だけ。通常の保存点再開とは別で、新Run・新継続・新封印を作らない。

Windowsの演算差や時間切れでLinuxの全Committed状態・対照トレース・対象進捗と一致しなければ、監査器は拒否する。候補の20番目の系も含め、失敗をそのまま残す。全上部撤去や全経路の受入へ自動的に進める手順ではない。

関連検査とソース再生成は新規出力先で実行できる。

```powershell
.venv/Scripts/python.exe "$Study/run_checks.py" --evidence $Raw --output "$Raw/checks"
if ($LASTEXITCODE -ne 0) { throw 'focused checks failed; preserve individual logs' }
.venv/Scripts/python.exe "$Study/check_source_reproduction.py" `
  --evidence $Raw --output "$Raw/source-reproduction"
if ($LASTEXITCODE -ne 0) { throw 'isolated source reproduction differed' }
```

## 未達項目

全上部撤去、全撤去中の材料・接触履歴、全解放後の凍結、独立塑性3D全経路、両実装の全経路刻み感度、完了までの途中再開は未達。既存小ケースの独立側128/256も約0.9394で降伏ドリフト制御に失敗している。無投影の既存参照へ架空の投影や新構成則を追加しない。P3b/P4各行・構造/水/JOINT全体・P6全体は未完了。

湿潤JOINTの8192/64、未完了4096/128全周期、大規模Tet10/MC/MCC/斜面SRM・仮想仕事上端・大規模降雨/性能/実CAD・長時間Richards等は、本手順のクラウド実行対象外。優先2の既存不合格値と全保護記録はREADMEのとおり保持する。追加クレジット購入・PR #13のマージは行わない。
