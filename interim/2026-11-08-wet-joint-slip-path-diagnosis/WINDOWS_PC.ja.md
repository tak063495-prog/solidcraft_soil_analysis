# Windows PC手順（今回未実行）

新先端のWindows CIとWindowsでの数値再現は未確認。過去CIを転用しない。本手順は既存54節点・12 Tet10湿潤JOINTだけを対象とし、低プレロードNCや範囲外の大規模計算を起動しない。完走しても全履歴の合格にはしない。

## 先端・環境・試験を記録

作業中の変更を退避できる既存チェックアウトで、PowerShellから行う。出力先は毎回未存在の名前にする。Rust 1.99.1、Python 3.13とロック済み参照依存を使う。既存の適合venvがあれば再作成しない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
if ($LASTEXITCODE -ne 0) { throw 'pull failed' }
git merge-base --is-ancestor d2b3dd4da9a1e92bbb6bc27b8d03acfc81f0903d HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
$D = 'interim/2026-11-08-wet-joint-slip-path-diagnosis'
$R = 'runs/NEW-windows-wet-slip'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
if (Test-Path $R) { throw 'use a new output directory' }
New-Item -ItemType Directory "$R/paths" | Out-Null
Start-Transcript -Path "$R/windows-transcript.log"
git show -s --format='%H %P %s' HEAD
Get-ComputerInfo | Select-Object WindowsProductName,WindowsVersion,OsBuildNumber
rustc -Vv
cargo -V
# 必要な場合だけ作成・インストールする
if (!(Test-Path '.venv/Scripts/python.exe')) { py -3.13 -m venv .venv }
.venv/Scripts/python.exe -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependencies failed' }
.venv/Scripts/python.exe -V
cargo build --locked --release *> "$R/build.log"
if ($LASTEXITCODE -ne 0) { throw 'build failed' }
cargo fmt --all --check *> "$R/fmt.log"
if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
cargo test --locked --release --test hydrostatic_volume --test joint_law --test joint_water *> "$R/rust-tests.log"
if ($LASTEXITCODE -ne 0) { throw 'Rust tests failed' }
foreach ($name in @('test_wet_joint_slip_path.py','test_hydrostatic_volume_audit.py','test_joint_water_reference.py','test_cut_volume_history.py')) {
  .venv/Scripts/python.exe -m unittest discover -s tests -p $name -v *> "$R/$name.log"
  if ($LASTEXITCODE -ne 0) { throw "Python failed: $name" }
}
Get-FileHash target/release/soil-fem-3d.exe -Algorithm SHA256
```

Rustは18件、Pythonは8+5+3+4=20件を期待するが、実際の実行数・除外・失敗・スキップ・終了コードをログから記録する。全体試験を行ったことにはしない。新先端CIの観測とは別である。

## 保存証拠だけのSHA照合と再計測

リポジトリに全新ZIPと旧 `wet-joint-finer-full-history` の19 ZIP等が揃っていることを確認する。以下は新規求解やWindows再開の試験ではない。保存済みLinux RunをPythonで比較する。

```powershell
.venv/Scripts/python.exe "$D/remeasure.py" --extract "$R/extracted" --output "$R/saved-remeasurement.json" *> "$R/remeasure.log"
if ($LASTEXITCODE -ne 0) { throw 'evidence mismatch; do not certify' }
```

全ZIP/エントリーSHA、工程端計測、共通受理点比較、元入力7再開Runの一致を再計測する。Linux実行版の封印をWindows実行版へ付け替えない。

## Windows実行版で新規求解

各実行は入力だけから生成する。独立参照は物理条件・既存Newton60・力1e-10・真の線形1e-10を保持する（生産基準の力1e-8より厳しい独立参照の元基準）。`--seconds` は外側の実時間上限であり、許容値・最小刻み・Newton上限ではない。PC性能によって上限終了し得る。

```powershell
.venv/Scripts/python.exe "$D/slip_study.py" native --solver target/release/soil-fem-3d.exe --output "$R/native-original" *> "$R/native-driver.log"
if ($LASTEXITCODE -ne 0) { throw 'native solve/restart/audit failed' }
# 元入力の全7工程・7個の厳密な途中保存点・全Run一致を新しいWindows封印で検査する。
.venv/Scripts/python.exe "$D/slip_study.py" bounded --steps 512 --local 64 --seconds 600 --output "$R/paths/g512-l64" *> "$R/g512-l64-driver.log"
.venv/Scripts/python.exe "$D/slip_study.py" bounded --steps 1024 --local 64 --seconds 900 --output "$R/paths/g1024-l64" *> "$R/g1024-l64-driver.log"
.venv/Scripts/python.exe "$D/slip_study.py" bounded --steps 512 --local 128 --seconds 900 --output "$R/paths/g512-l128" *> "$R/g512-l128-driver.log"
# クラウドで未完走の経路。まず1経路ずつ、予算を決めて実施する。
.venv/Scripts/python.exe "$D/slip_study.py" bounded --steps 4096 --local 128 --seconds 3600 --output "$R/paths/g4096-l128" *> "$R/g4096-l128-driver.log"
.venv/Scripts/python.exe "$D/slip_study.py" bounded --steps 8192 --local 64 --seconds 3600 --output "$R/paths/g8192-l64" *> "$R/g8192-l64-driver.log"
```

外側ドライバは未完了記録の保存に成功して終了コード0になり得る。**各経路の `terminal.json.complete` と `reason`、子プロセスの `returncode` を必ず読む。** `wall_clock_limit` は時間上限であり、Newton未収束や破壊ではない。Windowsの終了コードはLinuxの-15と一致するとは限らない。`status.json` と7工程/全受理列を確認し、結果ファイル欠落や途中工程だけを完走にしない。

独立参照のNPZは観測場であり再開できない。上限後は新しい別出力先で入力から再生成する。既存の途中状態を編集・再封印しない。7再開の対象は `native-original` であり、細分した独立経路のWindows再開成功を意味しない。

旧証拠は新規求解の後に読み、同じ比較器で計測する。

```powershell
.venv/Scripts/python.exe "$D/slip_study.py" saved --output "$R/saved" *> "$R/saved-driver.log"
if ($LASTEXITCODE -ne 0) { throw 'old evidence authentication failed' }
.venv/Scripts/python.exe "$D/slip_study.py" measure --output $R *> "$R/measure.log"
if ($LASTEXITCODE -ne 0) { throw 'measurement failed' }
.venv/Scripts/python.exe "$D/history_analysis.py" --evidence $R --output "$R/history-measurements.json" *> "$R/history-measure.log"
if ($LASTEXITCODE -ne 0) { throw 'common history audit failed' }
Get-Content "$R/measurements.json" -TotalCount 12
Stop-Transcript
Get-ChildItem $R -Recurse -File | Where-Object { $_.Name -ne 'files-sha256.csv' } | Get-FileHash -Algorithm SHA256 | Export-Csv "$R/files-sha256.csv" -NoTypeInformation
```

比較JSONは全7工程端の18場を測り、別JSONは保存した変位と接触履歴5場の共通受理点を比較する。未完了4096/128は保存接頭部の比較に限定し、欠けた工程を補間しない。8192/64が完走した場合でも、保存される比較は差の計測であり `passed=false` / `history_converged=false` を保持する。既存の全履歴受入基準がない状態で、新しい数値を合格閾値にしない。

低プレロードNCの0.9514000415802002・266増分と次の最小増分Newton60未収束、保護された他分野の結果は保持する。この手順からNC、大規模斜面/SRM、MCC空間再求解、降雨、実CAD、Native UI、GPU、完成配布へ対象を広げない。P3b/P4各行、構造・水・JOINT全体、P6全体は未完了のままである。
