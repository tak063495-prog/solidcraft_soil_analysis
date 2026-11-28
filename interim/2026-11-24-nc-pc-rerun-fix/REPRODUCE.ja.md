# PC 再実行・Linux 診断の再現

保存先の日付は依頼ラベル。実施日は 2026-11-12（日本時間）。リポジトリのルートから実行する。PR #13 はマージしない。出力先は毎回新しい名前にする。終了コード 0 は診断プロセスの正常終了であり、数値受入・全撤去・fold 確定を意味しない。NC-01 / P6 は未完のまま。

## 1. Windows のコードと環境

既存の b069ec8 の worktree 内で、作業中の変更を確認した上で実行する。別 worktree が同じブランチを保持している場合も、以下の detached checkout ならそのブランチを移動しない。以前の計算出力は保存する。

```powershell
git status --short
git fetch origin codex/plastic-spatial-convergence
git switch --detach origin/codex/plastic-spatial-convergence
git rev-parse HEAD
$Py = '.venv/Scripts/python.exe'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
& $Py -c "import sys,numpy,scipy; print(sys.version); print(numpy.__version__,scipy.__version__)"
```

固定依存は既存の `requirements-reference.lock.txt` と既存の Rust toolchain / Cargo.lock を使う。

## 2. 完了済み trace を再利用し、fold_path と analyze だけを再実行

`$TraceRun` は前回実行した `$TaskRun` の実際のフォルダーに変更する。`$FixRun` は未存在のフォルダーを指定する。`history/run.json`、`accepted.jsonl`、checkpoint を編集・再封印しない。前回の失敗出力も消さない。

```powershell
$TraceRun = 'runs/NEW-nc-fold-PC' # 前回の実際の保存先に変更する
$FixRun = 'runs/NEW-nc-pc-rerun-fix'
if (Test-Path $FixRun) { throw '新しい出力先を指定してください' }
New-Item -ItemType Directory -Path $FixRun -ErrorAction Stop | Out-Null
foreach ($denom in @(128,256)) {
  $p = "$TraceRun/series/$denom"
  if (-not (Test-Path "$p/history/mesh.json")) {
    # 元 mesh のコピーだけ。run / accepted / checkpoint は変更しない。
    Copy-Item "$p/mesh.json" "$p/history/mesh.json" -ErrorAction Stop
  }
  & $Py tools/nc_fixed_inputs.py --history "$p/history" --output "$FixRun/fixed-input-$denom.json"
  if ($LASTEXITCODE -ne 0) { throw "固定入力不一致: $denom" }
}
& $Py tools/prepare_nc_fold_observer.py --output "$FixRun/observer"
if ($LASTEXITCODE -ne 0) { throw '材料観測器のビルド失敗' }
$Probe = (Resolve-Path "$FixRun/observer/target/release/examples/material_diagnostic.exe").Path
foreach ($denom in @(128,256)) {
  $History = "$TraceRun/series/$denom/history"
  foreach ($scale in @('coarse','fine')) {
    $ds = if ($scale -eq 'coarse') { '0.00002' } else { '0.00001' }
    $count = if ($scale -eq 'coarse') { 128 } else { 256 }
    $PathOut = "$FixRun/native$denom-$scale"
    & $Py tools/nc_fold_path.py --map native --fresh-native $History --binary $Probe --step-size $ds --count $count --budget-seconds 1800 --output $PathOut
    if ($LASTEXITCODE -ne 0) { throw "fold_path 異常終了: $denom $scale （保存出力を保持）" }
    & $Py tools/analyze_nc_fold_path.py --path $PathOut --output "$PathOut-audit"
    if ($LASTEXITCODE -ne 0) { throw "analyze 異常終了: $denom $scale （保存出力を保持）" }
  }
}
```

`fixed-input-*.json` は、生の差分と正規化後の差分を保存する。今回のスキーマ差なら `raw_equal=false`、`fixed_inputs_match=true`、`canonical_differences=[]` になる。未知フィールドや実際の材料・荷重・制御変更は拒否する。履歴そのものは正規化した値で上書きしない。既存 trace の `complete=false` / `numerical_acceptance=false` も維持する。

この再実行には `nc_history_trace` の再ビルド・再求解は不要。既存履歴から read-only の材料観測器へ渡して局所経路を計算する。旧 Linux / Windows checkpoint を別バイナリで本番再開する手順ではない。局所経路の `report.json`、`fold-results.json` と各残差・荷重・場の監査を読み、点数不足・不受理・時間制限なら未確定を残す。直接経路と新しい診断 restart の比較は、同じ新しいソース・バイナリ・ds・入力で行う。

## 3. 線形残差を追加観測する場合だけ、新規前史を計算

既存 trace の再利用にこの節は不要。元行列・右辺・解を残して原因を診断するための追加計算である。既定求解器を変更せず、隔離コピーの observer が既存 `dense_lu` で拒否された線形系だけを記録する。異なるバイナリの checkpoint は再開しない。反復改良は保存行列に対するオフライン計算で、本番の受理判定には使わない。

```powershell
& $Py tools/prepare_nc_fold_series.py --output "$FixRun/linear-series"
if ($LASTEXITCODE -ne 0) { throw '入力生成失敗' }
& $Py tools/prepare_nc_linear_observer.py --output "$FixRun/linear-observer"
if ($LASTEXITCODE -ne 0) { throw '線形観測器ビルド失敗' }
$LinearTrace = (Resolve-Path "$FixRun/linear-observer/target/release/examples/nc_history_trace.exe").Path
$DenseProbe = (Resolve-Path "$FixRun/linear-observer/target/release/examples/nc_dense_lu_diagnostic.exe").Path
foreach ($denom in @(128,256)) {
  $p = "$FixRun/linear-series/$denom"
  $Systems = Join-Path (Resolve-Path $FixRun).Path "linear-$denom.jsonl"
  if (Test-Path $Systems) { throw '線形系の出力先が既に存在する' }
  try {
    $env:NC_LINEAR_SYSTEMS = $Systems
    & $LinearTrace direct "$p/mesh.json" "$p/ground.json" "$p/stages.json" "$FixRun/linear-native-$denom"
    if ($LASTEXITCODE -ne 0) { throw "線形観測の異常終了: $denom" }
  } finally { Remove-Item Env:NC_LINEAR_SYSTEMS -ErrorAction SilentlyContinue }
  & $Py tools/diagnose_nc_dense_lu.py --systems $Systems --binary $DenseProbe --refinements 8 --output "$FixRun/linear-audit-$denom"
  if ($LASTEXITCODE -ne 0) { throw "線形監査の異常終了: $denom" }
}
```

`diagnosis.json` には元の native 残差、同じ LU の再生、80 桁 Decimal による残差、後退誤差、行列積の丸め誤差尺度、条件数、各反復改良の解と残差を残す。許容値は `1e-10` のまま。最小値だけを示して失敗した試行を隠さない。反復改良で残差が下がっても、その解による非線形経路・全撤去・fold を受理したとは扱わない。

## 4. Linux の厳密な旧証拠再計測と全件検証

以下は今回の managed Linux 環境の例。固定依存の仮想環境と固定 Rust toolchain を有効にして実行する。

```bash
source /workspace/soil-environment/activate.sh
export RUST_TEST_THREADS=1
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets
cargo test --locked --release --doc
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence; do
  python tools/run_ci_python_tests.py --group "$group" --output "work/NEW-pc-rerun-validation/$group"
done
```

全群の発見 ID の一致、全 ID の欠落・重複なしの分割、failure / error / skip / ignored / filtered の件数を記録する。群が失敗しても残りの群を実行して記録する。失敗した実行を成功として数えない。

旧NCの2群を単独で再現する場合：

```bash
python interim/2026-11-10-nc-mean-stress-release/remeasure.py --output work/NEW-strict-mean
python interim/2026-11-10-nc-spatial-release/remeasure.py --output work/NEW-strict-spatial
```

各出力先は失敗しても保持する。Linux では `--allow-platform-float-differences` を指定しない。元アーカイブ、元 measurements、旧 hash / snapshot、許容値、既存拒否試験は変更しない。

全量再計測前に mean の最小例を診断するには、次を実行できる。`024835c` の元 source を Git から読み、元 SHA と照合する。shallow checkout 等で同 commit がない場合は、すでに SHA 検証して展開した元ファイルを明示する。

```bash
python tools/diagnose_nc_legacy_floats.py --output work/NEW-legacy-float-diagnosis
# 代わりに元展開ファイルを使う場合：
python tools/diagnose_nc_legacy_floats.py --fem-source work/NEW-strict-mean/evidence/sources/tools/fem_reference.py --output work/NEW-legacy-float-diagnosis-from-snapshot
```

この診断が終了コード 0 でも、元の厳密比較が合格した意味ではない。`diagnosis.json` の各差を確認する。旧 mean アーカイブは約 5.49 GB、spatial は約 0.92 GB。全件試験には全アーカイブの取得・検証・展開の時間と容量が必要となる。
