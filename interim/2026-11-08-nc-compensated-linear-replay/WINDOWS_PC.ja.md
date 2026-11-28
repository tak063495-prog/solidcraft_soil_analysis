# Windows PC：保存証拠の監査と限定線形再現

この先端のWindows CIは**未確認**。Linux Rust23件・Python17件の限定試験成功をWindows成功と扱わない。数値結果は `passed=false`、`history_converged=false`、0.9514000415802002・266増分、新規受理0件。最小1増分のNewton60失敗を破壊と呼ばない。

## ブランチと環境

既存出力を上書きしない。リポジトリルートからPowerShellで実行する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw 'branch update failed' }
git merge-base --is-ancestor eaeb160b68a1ec982634f45def65fc47bb279b45 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
git show -s --format='%H %P' eaeb160b68a1ec982634f45def65fc47bb279b45
# 親: 840de37cfc8b5ec7bab1fdf3ffed691c32a209ad
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$Study = 'interim/2026-11-08-nc-compensated-linear-replay'
```

Python仮想環境は既存 `requirements-reference.lock.txt` の固定依存を使う。再ビルド時だけ `rust-toolchain.toml` 指定のRust1.99.1が必要。以下では `.venv/Scripts/python.exe` を使う。OS、先端SHA、Python/Rust/依存バージョン、終了コードとログを新規出力先へ記録する。

## 推奨：求解せず、保存された証拠を検査

ZIPはPython標準zipfileで展開できるLZMA形式。新規名を使う。

```powershell
$Evidence = 'work/NEW-nc-linear-evidence'
$Output = 'work/NEW-nc-linear-audit'
if ((Test-Path $Evidence) -or (Test-Path $Output)) { throw 'choose fresh output paths' }
$Manifest = Get-Content "$Study/archive-sha256.json" -Raw | ConvertFrom-Json
$Actual = (Get-FileHash "$Study/evidence.zip" -Algorithm SHA256).Hash.ToLowerInvariant()
if ($Actual -ne $Manifest.archive_sha256) { throw 'archive checksum mismatch' }
.venv/Scripts/python.exe -m zipfile -e "$Study/evidence.zip" $Evidence
if ($LASTEXITCODE -ne 0) { throw 'extract failed' }
.venv/Scripts/python.exe "$Study/remeasure.py" --evidence $Evidence --output $Output
if ($LASTEXITCODE -ne 0) { throw 'evidence audit failed; retain failure' }
Get-Content "$Output/remeasurement.json"
```

全89エントリーのSHA256、19+60線形系の80桁監査、61 Newton場の保存応力からの力再組立、接触/撤去済み材料不変を検査する。Rust実行版・材料写像・NC本体は実行しない。Linuxでは4 JSONが完全一致した。WindowsのNumPy/BLASの差で完全一致しなければfalseを保存し、値を編集して一致させない。例外で拒否した場合も成功へ変更しない。`passed=false` と `history_converged=false` は維持する。

## 任意：隔離ソースで限定試験を再現

新規出力先に指定先端のソースを生成する。通常LU対照、19保存系の独立線形ゲート、直前の受理1増分対照、停止点直後の最小1増分、そこで生じる60線形系の監査だけを行う。全NCの初期状態からの実行や、正式な保存点再開は行わない。

```powershell
$Replay = 'work/NEW-nc-linear-replay'
if (Test-Path $Replay) { throw 'choose a fresh replay directory' }
.venv/Scripts/python.exe "$Study/reproduce.py" --output $Replay
if ($LASTEXITCODE -ne 0) { throw 'bounded reproduction differed or failed; preserve actual logs' }
Get-Content "$Replay/bounded-result.json"
```

候補非線形試験は、19/19保存系が独立ゲートを満たした場合だけ実行される。元の封印保存点・入力を改変せず、元実行版の識別子・型付きハッシュを検査する観測フックを用いる。異なるWindows実行版へLinux保存点の封印を付け替えない。通常の `resume-stages` の実行版検査を緩める手順はない。

Windowsの演算差で対照トレース、Committed状態、反復数が一致しなければ再現器は拒否する。実際の `diagnostic.json` とトレースを残し、Linux値を埋めて合わせない。診断プロセス終了コード0は収束・物理合格ではない。

関連ソフトウェア試験は次で行える。各コマンド後に `$LASTEXITCODE` を検査する。

```powershell
cargo fmt --all --check
cargo test --locked --release --lib linalg::tests
cargo test --locked --release --test free_contact
cargo test --locked --release --test mcc_midpoint
.venv/Scripts/python.exe "$Study/run_checks.py" python --output work/NEW-nc-linear-python-checks `
  --probe "$Replay/candidate/target/release/examples/nc_stop_probe.exe" `
  --prior "$Replay/evidence/prior" --traces "$Replay/evidence/diagnostic"
# reproduce.py が作成した隔離候補での追加3試験
cargo test --locked --release --manifest-path "$Replay/candidate/Cargo.toml" --test nc_linear_experiment
```

## 継続の条件と未達

候補は線形残差を改善したが、最小刻みの次増分は力残差約1.44e-7で元の1e-8を満たさない。今回の限定結果からNC全撤去を再試行する根拠は得られていない。力・真の線形・材料許容値、Newton上限、最小刻みを変更して通さない。

全撤去1.0、全撤去履歴、残存域λ=0、全解放と凍結、独立塑性3D全場、両実装の全経路刻み感度、完了までの別プロセス途中再開一致を揃えるまで受入しない。100 kPaプレロード合格を転用しない。P3b/P4各行、構造・水・JOINT全体、P6全体は未完。保護記録はREADME記載のまま保持し、指定大規模計算をクラウドで実行しない。PR #13をマージせず、追加クレジットを購入しない。
