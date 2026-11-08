# Windows：保存証拠の監査と小ケース診断の任意再実行

**全撤去は未達、passed=false。新先端のWindows CI・Windows実行は未確認。** この小ケースのクラウド診断にPC移行は必須ではない。以下は再現確認の手順であり、全解放や独立参照の受入手順ではない。

PowerShellでリポジトリルートから実行する。既存出力を上書きしない。既存 `requirements-reference.lock.txt` のPython環境、`rust-toolchain.toml` のRust 1.99.1を使用する。

```powershell
git switch codex/plastic-spatial-convergence
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw 'branch update failed' }
git merge-base --is-ancestor 06bff71dcc54f7060932aff57194f04ffe0593f8 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong branch history' }
git show -s --format='%H %P' 06bff71dcc54f7060932aff57194f04ffe0593f8
# 親: 6429e8366ad7b02ab546cc9a9f7ff321b693962b
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$Study = 'interim/2026-11-09-nc-turning-path-diagnosis'
```

## 推奨：求解せず再計測

```powershell
$Audit = 'work/NEW-nc-turning-audit'
if (Test-Path $Audit) { throw 'choose a fresh directory' }
.venv/Scripts/python.exe "$Study/remeasure.py" --output $Audit
if ($LASTEXITCODE -ne 0) { throw 'archive or measurements differed; retain actual output' }
Get-Content "$Audit/remeasurement.json"
```

両アーカイブと全エントリー、69ソース、元証拠の選択エントリーのSHA256を検査する。材料写像・Rustソルバーは起動しない。保存応力から釣合いを別再組立し、全128拡大線形系の真の残差を再計測する。新規展開したファイルだけへ書く。WindowsのNumPy/BLAS等の差で完全一致しなければ実際の差を保存し、Linux値を埋めない。

`all_matches=true` は保存計測の再現であり、数値受入ではない。`passed=false` と `history_converged=false` を維持する。

## 任意：隔離ソースから診断を新規実行

```powershell
$Raw = 'work/NEW-nc-turning-run'
if (Test-Path $Raw) { throw 'choose a fresh directory' }
.venv/Scripts/python.exe "$Study/capture.py" --output $Raw
if ($LASTEXITCODE -ne 0) { throw 'diagnostic failed; retain real results and logs' }
Get-Content "$Raw/audit/measurements.json"
Get-Content "$Raw/final-binary-sha256.json"
.venv/Scripts/python.exe "$Study/test_audit.py" `
  --evidence $Raw --output "$Raw/guard-tests" -v
if ($LASTEXITCODE -ne 0) { throw 'audit guard tests failed' }
```

元アーカイブ1,710エントリーを照合してから、指定06bff71dの数値ソースに観測関数だけを加えた隔離crateを新規生成する。固定分割材料診断3,240件と、係数の逆行を許す変位制御診断64点×2経路を実行する。本体src・通常保存点・元履歴・封印を編集せず、新Runや新継続を作らない。Windowsの実行版SHAと実測出力をそのまま保存する。LinuxのバイナリSHAへの付け替えや通常再開はしない。

係数が折り返す診断経路は、指定された単調撤去ではない。最後まで64点が求まっても撤去1.0・破壊・全履歴収束・独立塑性全場・完了再開の合格にはしない。途中エラーや時間切れは、その理由と到達点を記録する。

## Linuxで実行した範囲に相当する関連検査

```powershell
cargo fmt --all --check
if ($LASTEXITCODE -ne 0) { throw 'fmt failed' }
cargo test --locked --release --test mcc_midpoint --test free_contact
if ($LASTEXITCODE -ne 0) { throw 'selected Rust tests failed' }
cargo clippy --locked --release --manifest-path "$Raw/isolated/Cargo.toml" `
  --target-dir "$Raw/target" --example material_diagnostic -- -D warnings
if ($LASTEXITCODE -ne 0) { throw 'isolated clippy failed' }
```

Linux実績は関連Rust12件、Python監査7件。全Rust/PythonスイートやWindows CIを実行したという意味ではない。新しい物理仕様・形状・参照値・合格値は追加しない。追加クレジット購入、PR #13のマージ、大規模求解をこの手順に含めない。
