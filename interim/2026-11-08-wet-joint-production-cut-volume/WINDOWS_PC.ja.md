# Windows PCでの照合・再計測・限定再求解

このフォルダの `passed=false` / `history_converged=false` は終了コード0で変更しない。新先端のWindows CI成功とWindowsでの数値再現は未確認。PR #13をマージしない。追加クレジットを買わない。以下はリポジトリ直下のPowerShellで実行する。

## ブランチと依存関係

```powershell
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 840de37cfc8b5ec7bab1fdf3ffed691c32a209ad HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
git status --short
git log -1 --format='%H %P'
git ls-remote origin refs/heads/codex/plastic-spatial-convergence refs/pull/13/head
```

ローカル先端とGitHubの同一性を確認する。指定起点の親は `09ff926d82e574561393be2ac877cf9f5627d750`。ユーザーPCで変更がある場合は上書きしない。

既存の固定依存関係と `rust-toolchain.toml`（Rust 1.99.1）を使用する。

```powershell
py -3.13 -m venv .venv
.venv/Scripts/python.exe -m pip install --require-hashes -r requirements-verification-win-cp313.lock.txt
$StudyRoot = 'interim/2026-11-08-wet-joint-production-cut-volume'
$RunRoot = 'runs/NEW-wet-production-cut'
```

毎回、まだ存在しない出力名を選ぶ。

## 求解なしで生証拠を検査

```powershell
Get-FileHash "$StudyRoot/evidence.zip" -Algorithm SHA256
.venv/Scripts/python.exe "$StudyRoot/verify_evidence.py" --output "$RunRoot-readback"
if ($LASTEXITCODE -ne 0) { throw 'evidence verification failed' }
```

`archive-sha256.json` のZIP全体と全エントリーを照合し、新規領域へ展開する。既存1024/32・2048/32のZIPと、2048/64・4096/64の19分割ZIPを元のリポジトリの場所に全て揃える。独立参照を新しく求解せず、全18場・保持荷重・保存点の力監査・全計測を再現する。成功はソフトウェア再計測の一致であって全履歴受入ではない。

## ビルドと既存回帰試験

```powershell
cargo fmt --all --check
cargo clippy --locked --release --all-targets -- -D warnings
cargo test --locked --release
cargo build --locked --release
cargo build --locked --release --example hydrostatic_volume_probe
.venv/Scripts/python.exe "$StudyRoot/run_python_tests.py" --output "$RunRoot-python-tests.json"
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v
.venv/Scripts/python.exe "$StudyRoot/check_integration.py" --probe target/release/examples/hydrostatic_volume_probe.exe --output "$RunRoot-integration"
```

各コマンドの終了コードを確認する。Pythonランナーは低プレロードNCを実際に起動するクラス6件と別メソッド1件を明示的に除外し、一覧・件数を保存する。これらを「成功」「スキップ成功」と数えない。入力検査・既存保存場の読出しだけのNC試験は実行する。新しい改善根拠がないため、この作業のために低プレロードNCを再試行しない。

## 同じ小ケースだけを再求解

```powershell
.venv/Scripts/python.exe "$StudyRoot/study.py" --solver target/release/soil-fem-3d.exe --output "$RunRoot-solve"
if ($LASTEXITCODE -ne 0) { throw 'bounded wet JOINT study failed; inspect logs' }
```

元入力と256/512分割・接触最大刻み5e-6/2.5e-6の4生産経路、256分割の全7工程それぞれの途中再開を実行する。形状・構成則・物理荷重・水位・力/線形許容値は変更しない。生産接触分割は適応的で、独立参照の固定32/64とは別。

Linux保存点はLinux実行版のSHA256に封印されている。Windows実行版へ封印を書き換えたり、Linux保存点をそのままWindowsで再開したりしない。上のコマンドはWindowsの同じ実行版で新規保存点と連続実行を作り、そのRun全体の一致を確認する。

生証拠の数値と別OSでのビット一致は仮定しない。全7工程の完了、全受理点の力1e-8・真の線形1e-10、元GP表示、物理水圧2.5 kPa、独立力再組立、保持荷重、全18場、全体/接触刻み感度、途中再開一致を確認する。今回の独立参照器を一般塑性・部材・Darcy/Richards・水位変更・掘削・BIRTHへ拡張しない。全履歴の既存基準がない項目は不合格/未検証のまま残す。

P3b/P4各行、構造・水・JOINT全体、P6全体は未完了。Griffiths1200/4800/19200 Tet10、MCC274〜14909、斜面SRM本計算、82944降雨、105456性能、実CAD、長時間Richards等をこの限定手順で起動しない。
