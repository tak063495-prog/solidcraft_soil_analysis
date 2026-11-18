# PC再実行手順（この作業ではPC計算を実行しない）

Windows / GPD WIN Max 2 / 利用可能メモリ約48GB / `D:\Work` / PowerShell 7。
対象PR #13、`codex/plastic-spatial-convergence`。マージ・PCからのpushは行わない。
既存のPID 18836・旧worktree・旧失敗出力を操作しない。以下は全て新しいworktree/一意の出力。

実行予算は各子7500秒以下、原則8GiB以下（証拠/CLI4GiB）、1スレッド。
`pcfail_run.py`はWindows Job Objectへ停止状態の新規子を所属させてから再開する。
既存PIDを探して停止しない。Job総コミットメモリを制限し、旧ログのpeak working setとは別の量を保存する。
タイムアウトではそのJobだけ終了。Job割当失敗は未実行失敗とし、無制限実行に降格しない。
Windows Job Object経路は本作業のLinuxでは実行未確認、push後Windows CIも未確定。
不具合時は出力を保存して停止する。全体上限40GB/1経路1日まで拡大して続ける手順ではない。

## 0. sparse worktree、環境、ビルド

cwd：`D:\Work`。初回取得もblob:none / depth1 / single-branch。
既存のNEW-*リポジトリの設定や実行中作業ツリーを変えず、別のオブジェクト格納先を作る。
`$Commit`はPC再現に必要なコード・科学証拠の固定SHA。後続の台帳/検証記録追加とは分け、枝が動いても同じSHAを用いる。

```powershell
Set-Location D:\Work
$ErrorActionPreference='Stop'
$Stamp=(Get-Date -Format 'yyyyMMdd-HHmmss')+'-'+[guid]::NewGuid().ToString('N').Substring(0,8)
$Store="D:\Work\PCFAIL-objects-$Stamp"
$Tree="D:\Work\PCFAIL-$Stamp"
$Branch='codex/plastic-spatial-convergence'
$Commit='ce8ce8ddba71fba7fe4015696268669b6df72886'
git clone --filter=blob:none --depth=1 --single-branch --branch $Branch --no-checkout https://github.com/tak063495-prog/solidcraft_soil_analysis.git $Store
if ($LASTEXITCODE -ne 0) { throw '限定clone失敗' }
git -C $Store fetch --filter=blob:none --depth=1 --no-tags origin $Commit
if ($LASTEXITCODE -ne 0) { throw '固定commit取得失敗' }
git -C $Store worktree add --detach --no-checkout $Tree $Commit
if ($LASTEXITCODE -ne 0) { throw '別worktree作成失敗' }
git -C $Tree config core.autocrlf false
git -C $Tree sparse-checkout init --no-cone
$Patterns=@('/*','!/interim/**/*.zip','!/interim/**/*.zip.part*','!/interim/**/*.tar',
 '!/interim/**/*.tar.gz','!/interim/**/*.tar.part*','!/interim/**/*.tar.gz.part*',
 '!/interim/**/*.tgz','!/interim/**/*.gz','!/interim/**/*.7z','!/interim/**/*.npz','!/interim/**/*.part*')
$Patterns | git -C $Tree sparse-checkout set --no-cone --stdin
if ($LASTEXITCODE -ne 0) { throw 'sparse設定失敗' }
git -C $Tree read-tree --reset -u HEAD
if ($LASTEXITCODE -ne 0) { throw 'sparseソース配置失敗' }
Set-Location $Tree
if ((git rev-parse HEAD).Trim() -ne $Commit) { throw '固定SHA不一致' }
$Diag='interim/2027-01-08-pc-failure-diagnosis'
$Py="$Tree\.venv\Scripts\python.exe"
$Out="$Tree\runs\pcfail-$Stamp"
New-Item -ItemType Directory $Out | Out-Null
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
py -3.13 -m venv .venv
if ($LASTEXITCODE -ne 0) { throw 'Python3.13 venv失敗' }
& $Py -m pip install --disable-pip-version-check -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の導入失敗' }
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
if ($LASTEXITCODE -ne 0) { throw 'Rust1.99.1導入失敗' }
rustc --version --verbose
& $Py -c "import numpy as n; print('longdouble',n.finfo(n.longdouble),'float64',n.finfo(n.float64))" | Out-File "$Out\precision.txt" -Encoding utf8
& $Py "$Diag/pcfail_run.py" --output "$Out/build-process" --timeout 7500 --memory-gib 8 -- cargo build --locked --release --all-targets
if ($LASTEXITCODE -ne 0) { throw 'build未完了。新規出力を保持' }
```

所要：取得は回線次第、ビルド10–60分/概算2–8GiB（PC未実測）。
source-ishの規模は依頼上312–461MB。Cargo/venv/targetは別。
29GBのフルcheckout・全履歴取得・全アーカイブ展開は不要。
見る値：HEAD、Rust1.99.1、Python3.13、longdouble精度、build-processの終了/時間/Jobメモリ。
不足依存・固定blob取得失敗はincomplete/未達として残す。

## 1. 証拠partだけを取得・検証（b01-evidence / a01-evidence）

cwd：`$Tree`。先行条件：手順0のソース配置とvenv。
取得前に欠落を記録する場合は次を実行する。終了2がincompleteであり、合格/skipではない。
独立した別項目は手動で続けてよいが、この証拠の再集計へは進まない。

```powershell
& $Py "$Diag/pcfail_evidence.py" --profile b01 --output "$Out/b01-before-fetch"
# 終了2: evidence-status.jsonの全missing partを保持。終了1:破損等を保持して原因確認。
```

次は**必要な4群12part、計207,550,464 bytesだけ**を固定commitからblob選択取得する。
除外パターンは維持し、ソースと証拠を混ぜた全面checkoutはしない。
各partのサイズ/SHAを先に確認し、既存破損partは上書きしない。
取得失敗やSHA不一致で止め、既存ファイルを削除する手順はない。

```powershell
foreach ($EvidenceProfile in @('b01','b01-cli','a01','a01-old')) {
  & $Py "$Diag/pcfail_run.py" --output "$Out/fetch-$EvidenceProfile-process" --timeout 1800 --memory-gib 4 -- `
    $Py "$Diag/pcfail_evidence.py" --profile $EvidenceProfile --fetch-commit $Commit --output "$Out/fetch-$EvidenceProfile"
  if ($LASTEXITCODE -ne 0) { Write-Host "$EvidenceProfile の証拠は未達。記録を保持し、次の独立群を確認する" }
}
& $Py "$Diag/pcfail_run.py" --output "$Out/b01-evidence-process" --timeout 1800 --memory-gib 4 -- `
  $Py "$Diag/pcfail_evidence.py" --profile b01 --output "$Out/b01-evidence" --extract "$Out/b01-unsealed"
& $Py "$Diag/pcfail_run.py" --output "$Out/a01-verify-process" --timeout 1800 --memory-gib 8 -- `
  $Py "$Diag/pcfail_evidence.py" --profile a01 --output "$Out/a01-verify"
& $Py "$Diag/pcfail_run.py" --output "$Out/a01-old-process" --timeout 3600 --memory-gib 8 -- `
  $Py "$Diag/pcfail_evidence.py" --profile a01-old --output "$Out/a01-old" --extract "$Out/a01-old-subset"
```

群別：b01=2part/20,380,981B、b01-cli=4part/51,221,588B、a01=1part/21,860,534B、
a01-old=5part/114,087,361B。各1–5分、概算1GiB未満、回線時間は別。
検証は一時結合をディスクへ流し、全entryをSHA検証。b01展開58,390,662B、a01-oldは旧手順と同じ
selected/interacting-stacked-attemptの非checkpoint・非restarts部分だけ。全27.66GBを展開しない。
見る値：各evidence-status.jsonのstatus=verified、all_bytes_verified、part/全entry SHA、
passed=false、numerical_acceptance=false。証拠検証成功は元の数値受入を意味しない。

b01再集計はb01-evidence-processが終了0の場合だけ：

```powershell
& $Py "$Diag/pcfail_run.py" --output "$Out/b01-remeasure-process" --timeout 1800 --memory-gib 4 -- `
  $Py interim/2026-12-25-b01-curved-tri6-cause-diagnosis/b01cause_report.py `
  --formal "$Out/b01-unsealed/b01cause-formal" --output "$Out/b01-remeasured.json"
```

## 2. b01-cli（旧21拒否＋4受理対照）

cwd：`$Tree`。先行条件：build成功、b01-cliの4partがverified。
欠落時は入口がincompleteを記録し、solverを起動しない。

```powershell
& $Py "$Diag/pcfail_run.py" --output "$Out/b01-cli-process" --timeout 1800 --memory-gib 4 -- `
  $Py "$Diag/pcfail_cli.py" --solver "$Tree/target/release/soil-fem-3d.exe" --output "$Out/b01-cli"
```

所要1–5分/概算2GiB未満。出力：b01-cli/evidence-status.json、controls/cli-controls-observations.json。
見る値：refusals=21、accepted_controls=4、all_original_contracts_confirmed、入力SHA・終了コード・診断文字・
受理Runの旧比較値。曲面受入や真のcheckpoint再開の証拠へ読み替えない。

## 3. b02-n4、条件付きn8、さらに条件付きn16

cwd：`$Tree`。先行条件：build成功。旧driverでなくpcfail_b02.pyを使う（旧ファイルは失敗記録とともに保持）。
既存公開API観測器を別フォルダでビルドする。

```powershell
& $Py "$Diag/pcfail_run.py" --output "$Out/b02-observer-build-process" --timeout 7500 --memory-gib 8 -- `
  $Py interim/2026-12-01-tri6-offset-composite-b02/build_observer.py "$Tree/work/pcfail-b02-observer-$Stamp"
if ($LASTEXITCODE -ne 0) { throw '観測器build未完了' }
& $Py "$Diag/pcfail_run.py" --output "$Out/b02-n4-process" --timeout 7500 --memory-gib 8 -- `
  $Py "$Diag/pcfail_b02.py" --output "$Out/b02-n4" --observer "$Tree/target/release/b02-observer.exe" `
  --cli "$Tree/target/release/soil-fem-3d.exe" --steps 4
```

12/96 Tet10だけ。n4は10–60分/概算1–6GiB、7500s/8GiBを上限とする。PC実測ではない。
見る値：plan-results.jsonの2case、全restartsのexact_except_elapsed、native/plain/independentの終了0、
fresh_popen/invocation_id、Run全体と後続accepted列のexact（elapsedのみ除外）、元の各比較ゲート。
同じPIDが後で現れても失敗にしない。別起動の同一invocation_id・状態改変・履歴改変は失敗。

**次は自動チェーンへ登録しない。** n4の全工程・厳密再開検査が終了0で、全caseのcomparison.existing_field_gates_passed=true、
1子の実測時間/メモリから次も7500s/8GiB内と判断できた場合だけ、別コマンドとしてn8を選ぶ。
n8完了と同じ資源条件を確認してからn16。上限外の見込みなら規模上限のため確認不可・未達。

```powershell
# 条件を満たしたときだけ、まずn8（目安20–90分/1–6GiB）。
& $Py "$Diag/pcfail_run.py" --output "$Out/b02-n8-process" --timeout 7500 --memory-gib 8 -- `
  $Py "$Diag/pcfail_b02.py" --output "$Out/b02-n8" --observer "$Tree/target/release/b02-observer.exe" `
  --cli "$Tree/target/release/soil-fem-3d.exe" --steps 8
# n8が上の条件を満たした場合だけn16（目安30–120分/1–6GiB）。
& $Py "$Diag/pcfail_run.py" --output "$Out/b02-n16-process" --timeout 7500 --memory-gib 8 -- `
  $Py "$Diag/pcfail_b02.py" --output "$Out/b02-n16" --observer "$Tree/target/release/b02-observer.exe" `
  --cli "$Tree/target/release/soil-fem-3d.exe" --steps 16
```

既存b02cause_report.pyは4/8/16が同じrootにある契約のため、上記の別rootへそのまま呼ばない。
個別plan-resultsと生出力を保持する。今回の課題はPCのn4を止めたassertの修正で、
別rootを自動連結して未実行系列を実行済み扱いにする処理は加えていない。

## 4. b01-roof-n8（再失敗も結果として保存）

cwd：`$Tree`。先行条件：build成功。許容値・cutbackは旧入力のまま。
今回クラウドでもn8は未完了なので、これは修正済み経路の合格確認ではなくWindows再現。

```powershell
$Roof="$Out/roof-n8"
& $Py interim/2026-12-03-tri6-curved-b01/b01_cases.py $Roof --kind roof --n 8 --steps 4
if ($LASTEXITCODE -ne 0) { throw '入力生成失敗' }
& $Py "$Diag/pcfail_run.py" --output "$Out/roof-n8-process" --timeout 1800 --memory-gib 8 -- `
  "$Tree/target/release/soil-fem-3d.exe" solve-stages "$Roof/mesh.json" "$Roof/ground.json" "$Roof/stages.json" "$Roof/run.json"
```

目安1–5分/概算1GiB未満（旧PC求解実測4.2s/peak_ws113,528,832B、クラウドはREADME参照）。
見る値：run.complete、frames/attempts、last committed、全rejection reason、state unchanged。
この128平面Tri6＋12土Tet10以外へ進めない。n8の失敗を曲面ガード拒否・物理破壊へ変更しない。
観測用の別Rustビルドと行列解析はクラウドで保存済み。PCでは通常CLIだけでよい。

## 5. a04-matrix、条件付き32/64/128

cwd：`$Tree`。先行条件：build成功（a04_checkpoint_rehashを含む）。
独立参照の精度不足を許容値変更で回避しない。既存a04cause_driverをそのまま実行し、
failure.json・soil.log・独立部分出力を保持する。Windows binary64で再失敗する可能性がある。

```powershell
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-matrix-process" --timeout 7500 --memory-gib 8 -- `
  $Py interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis/a04cause_driver.py --output "$Out/a04-matrix"
# 原因観測だけを追加する場合。元の反復条件を変えない独立入力専用診断。
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-precision-process" --timeout 1800 --memory-gib 4 -- `
  $Py "$Diag/pcfail_hydraulic.py" --output "$Out/a04-precision" --steps 4 --precision host
```

本作業はportable高精度参照を実装していない。hostのnmantが52でmatrix未完了なら、32以降へ進まない。
将来参照を修正するときは、その修正を含む検証済み固定SHAを別記し、今回の旧失敗出力を保持する。

matrix目安1–10分/概算1GiB未満（旧PC18.0s、peak_ws170,000,384Bで途中停止）。
精度観測は1分以内/概算1GiB未満。observation.jsonのexecution_complete=falseでも観測コマンド自体は
記録完了として終了0になりうる。これはmatrix成功条件に使わない。
見る値：driver-result.execution_complete、全native/独立exit_code、各comparison/元ゲート、
failure.json、hydraulic line search、longdoubleのmantissa、未完了/保存状態数。

**matrixの本体・独立・audit・convergenceがすべて終了0かつ元の限定比較を満たし、
実測から次が7500s/8GiB内の場合だけ**32を選ぶ。32の同条件を確認して64、64後に128。
現時点でWindows matrixは未達なので、下は条件付きコマンドだけで自動再実行対象ではない。

```powershell
# まず32のみ（目安5–20分/1–2GiB）。
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-32-process" --timeout 7500 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/study.py --output "$Out/a04-32" --steps 32
if ($LASTEXITCODE -ne 0) { throw 'a04-32未完了。次水準へ進まない' }
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-32-audit-process" --timeout 1800 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/audit.py --run "$Out/a04-32"
if ($LASTEXITCODE -ne 0) { throw 'a04-32の元ゲート監査未達。次水準へ進まない' }
# 32完走・元ゲート・予算を確認後だけ64（目安10–40分/1–2GiB）。
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-64-process" --timeout 7500 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/study.py --output "$Out/a04-64" --steps 64
if ($LASTEXITCODE -ne 0) { throw 'a04-64未完了。次水準へ進まない' }
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-64-audit-process" --timeout 1800 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/audit.py --run "$Out/a04-64"
if ($LASTEXITCODE -ne 0) { throw 'a04-64の元ゲート監査未達。次水準へ進まない' }
# 64の同条件確認後だけ128（目安20–80分/1–3GiB）。
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-128-process" --timeout 7500 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/study.py --output "$Out/a04-128" --steps 128
if ($LASTEXITCODE -ne 0) { throw 'a04-128未完了。次水準へ進まない' }
& $Py "$Diag/pcfail_run.py" --output "$Out/a04-128-audit-process" --timeout 1800 --memory-gib 8 -- `
  $Py interim/2026-11-29-richards-wet-joint-a04/audit.py --run "$Out/a04-128"
if ($LASTEXITCODE -ne 0) { throw 'a04-128の元ゲート監査未達。次水準へ進まない' }
```

これらは同じ12 Tet4だけ。見積りは保守的計画値でPC実測ではない。
3水準必須のa04cause_driver/convergence.pyへ単一の--stepsを渡さない。
各単独studyは元の本体/独立比較・checkpoint・auditを保ち、異なる刻み間の収束率はここでは算出しない。
行列完走前にdtを変えてline search失敗を消すことも、参照の厳しさを下げることも行わない。

## 6. 新しい軽量回帰試験と結果の読み方

この試験は今回の封印も読むため、基本sparseで除外される今回の4part（61,109,258B）だけを
追加で選択取得する。旧4群12partとは別。汎用fetch APIはmanifestのpartだけを固定blobから取得し、
sparse規則の全面再適用や旧失敗出力の削除を行わない。

```powershell
$NewEvidenceCode="import sys,json; from pathlib import Path; sys.path.insert(0,'$Diag'); from pcfail_evidence import fetch,verify; f=Path('$Diag').resolve(); m=json.loads((f/'archive-index.json').read_text(encoding='utf-8')); fetch(f,m,'$Commit'); r=verify(f,m); print(r); assert r['status']=='verified'"
& $Py "$Diag/pcfail_run.py" --output "$Out/new-evidence-process" --timeout 1800 --memory-gib 4 -- `
  $Py -c $NewEvidenceCode
if ($LASTEXITCODE -ne 0) { throw '新試験用の4partが未達。欠落を試験成功へ変換しない' }
& $Py "$Diag/pcfail_run.py" --output "$Out/new-tests-process" --timeout 1800 --memory-gib 4 -- `
  $Py -m unittest discover -s tests -p test_pc_failure_diagnosis.py -v
```

cwd `$Tree`、数秒/1GiB未満。PID再利用、厳密状態改変拒否、missing/corrupt、選択blob復元、
独立続行・timeout、観測用コピー、実CRLF末尾/OS非依存パスを検査する。
既存全Python試験は旧アーカイブを多数要求するため、このPC手順で全partを追加して29GBへ戻さない。
全件は本作業のクラウドvalidation.jsonを参照する。新試験成功を全件成功へ読み替えない。
全ソフトウェア試験の件数と各数値課題の未達を分ける。PCからpushせず、旧失敗出力を消さない。
