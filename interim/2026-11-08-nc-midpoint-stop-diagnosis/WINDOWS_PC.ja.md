# Windows PC：低プレロードNC停止点の診断・保存証拠の再計測

今回の新先端に対するWindows CIは**未確認**。以下はユーザーPC向け手順で、Windowsでの成功を記録したものではない。今回の結論は `passed=false`、`history_converged=false`。0.9514000415802002・266増分からの全撤去再開は行っていない。CPUの高速化・時間予算の延長だけで通るという根拠もない。

## 取得と環境

リポジトリルートで実行する。元の履歴・出力・保存点を上書きしない。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw 'branch update failed' }
git merge-base --is-ancestor 09ff926d82e574561393be2ac877cf9f5627d750 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong history' }
git show -s --format='%H %P' 09ff926d82e574561393be2ac877cf9f5627d750
# 親は c4de9fdafa4d7147f4f2bd8c8c877ca969b269b6
```

Python仮想環境へ既存 `requirements-reference.lock.txt` の固定依存を用意する。Rustの診断実行版をビルドする場合は既存 `rust-toolchain.toml` の1.99.1を用意する。先端SHA、OS、Python/Rustと依存バージョン、各コマンドの終了コードとログを別々に保存する。

## 推奨：求解せず保存証拠を再計測する

未存在の出力名を使う。ZIPはLZMA方式で、Python標準zipfileで展開できる。

```powershell
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$Study = 'interim/2026-11-08-nc-midpoint-stop-diagnosis'
$Evidence = 'work/NEW-nc-stop-evidence'
$Output = 'work/NEW-nc-stop-remeasured'
$Manifest = Get-Content "$Study/archive-sha256.json" -Raw | ConvertFrom-Json
$Actual = (Get-FileHash "$Study/evidence.zip" -Algorithm SHA256).Hash.ToLowerInvariant()
if ($Actual -ne $Manifest.archive_sha256) { throw 'archive checksum mismatch' }
if ((Test-Path $Evidence) -or (Test-Path $Output)) { throw 'choose fresh output paths' }
.venv/Scripts/python.exe -m zipfile -e "$Study/evidence.zip" $Evidence
if ($LASTEXITCODE -ne 0) { throw 'extract failed' }
.venv/Scripts/python.exe "$Study/remeasure.py" --evidence $Evidence --output $Output
if ($LASTEXITCODE -ne 0) { throw 'saved-evidence audit failed; no acceptance' }
Get-Content "$Output/remeasurement.json"
```

この処理は全エントリーSHAを検査し、保存された20 Newton場・48線探索候補と材料応答を読む。保存点を変更・再封印せず、CLI・材料写像・独立全経路を実行しない。残存域λ=0、保存点の力残差、接触履歴・削除材料の不変、3240材料応答の一致、接線のSVD等を再計測する。

`solver_executed=false`、`material_map_executed=false`、`passed=false` が期待される。`exact_saved_json_matches` は保存計測との完全一致を個別に示す。WindowsのBLAS・long-double等の差で浮動小数値が異なる場合はfalseをそのまま保存し、値を編集して合わせない。既存の力/状態監査が例外で拒否した場合も成功に変えない。Linuxでの4 JSON完全一致は `remeasurement.json` に記録しており、Windows一致の保証ではない。

## 任意：隔離コピーで記録済み2増分だけを観測する

これは継続計算・全撤去再開ではない。停止機構を調べる場合だけ使い、同条件の本計算を時間だけ延ばして繰り返さない。生産ファイルは編集せず、元コミットの数値ソースとこのフォルダの観測コードを新規 `work/` へ準備する。

```powershell
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$Study = 'interim/2026-11-08-nc-midpoint-stop-diagnosis'
$Diagnostic = 'work/NEW-nc-recorded-trials'
if (Test-Path $Diagnostic) { throw 'choose a fresh diagnostic path' }
.venv/Scripts/python.exe "$Study/capture.py" --output $Diagnostic
if ($LASTEXITCODE -ne 0) { throw 'diagnostic failed; no acceptance' }
Get-Content "$Diagnostic/measurements/measurements.json"
.venv/Scripts/python.exe "$Study/test_diagnosis.py" `
  --probe "$Diagnostic/source/target/release/examples/nc_stop_probe.exe" `
  --prior "$Diagnostic/prior" --traces "$Diagnostic/diagnostic" -v
if ($LASTEXITCODE -ne 0) { throw 'diagnostic guard tests failed' }
```

出力は `prior-binding.json`、隔離 `source/`、`build.log`、`diagnostic/*/trace.jsonl`・`diagnostic.json`、`process.json`、`binary-sha256.json`、`measurements/`。元NCアーカイブをSHA検査して取り込み、記録された最後の受理1増分と最小刻み拒否1増分だけを評価する。改変された入力/履歴/元実行版識別子・観測窓外の進捗は拒否する。Windowsの演算差でLinuxの記録状態・停止値に完全一致しなければ、診断の比較は拒否/未一致となる。これを数値成功としない。

この観測フックは保存点の実行版封印を変更する再開APIではない。異なるWindows実行版へLinux保存点の封印を付け替えない。新しい継続Runや保存点は生成しない。終了コード0も物理合格ではない。

関連するソフトウェア検証は次で行える。今回Linuxで実行した範囲はfmt、隔離exampleのclippy、既存関連Python9件、新規ガード5件であり、全Rust/Python試験成功へ拡大して記録しない。

```powershell
cargo fmt --all --check
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_nc_midpoint_resume.py' -v
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_low_preload_nc_diagnosis.py' -v
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_nc_newton_trace.py' -v
```

各コマンド後に `$LASTEXITCODE` を確認する。`validation.json` はLinuxの実測であり、Windows CIの結果ではない。

## 今後の再開条件と未達

現状の診断だけでは全撤去へ通る改善を実証していないため、今回の手順は全撤去を続行しない。検証済み改善が別に得られた場合だけ、既存の正しい入力/実行版結合と保存点の契約を守って途中再開を検討する。単一のNumPy線形解が1e-10を通ることや、1増分の保存状態一致を全経路の改善根拠にしない。

対象物理条件は必ず5 kPa追加プレロード、pc=100 kPa、1/0.5/0 kN。100 kPaプレロードの合格を転用せず、力1e-8・真の線形1e-10等を緩めない。全上部撤去1.0、保持荷重ゼロ、撤去中の材料/接触全履歴、残存域λ=0、全解放後の凍結、入力だけからの独立塑性3D全場、両実装の全経路刻み感度、完了までの別プロセス途中再開一致が揃うまで合格としない。時間切れ・未収束・切戻し枯渇を破壊と呼ばない。

P3b/P4各行、構造・水・JOINT全体、P6全体は未完。保護記録は [README.ja.md](README.ja.md) のとおり元のまま保持する。新しい亀裂先端・接触BIRTH・高次板・弱結合仕様を作らず、指定された大規模計算をクラウドで実行しない。PR #13をマージせず、追加クレジットを購入しない。
