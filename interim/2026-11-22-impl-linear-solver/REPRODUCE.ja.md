# NC 96 Tet10 の再実行

元 PC の `fine_20261122/nc96` に保存した **初期入力** 3 件を使うと、2026-11-22 の実行と比較できる。過去の出力・チェックポイントを保存したまま、新規ディレクトリで各経路を初期状態から解く。PC の入力・行列の SHA は T9 に転記されていないため、クラウド生成入力との同一性は実ファイルで確認する。

## Windows での準備

PR #13 を取得したリポジトリのルートで実行する。Rust は rust-toolchain.toml 指定版、生成手順の Python は 3.13 と requirements-reference.lock.txt 指定のパッケージを使う。出力名は未使用のものを指定する。

```powershell
git fetch origin codex/plastic-spatial-convergence
if ($LASTEXITCODE) { throw 'fetch failed' }
git switch --detach origin/codex/plastic-spatial-convergence
if ($LASTEXITCODE) { throw 'checkout failed' }
git rev-parse HEAD
cargo build --locked --release --example nc_history_trace
if ($LASTEXITCODE) { throw 'build failed' }
$Trace = (Resolve-Path target/release/examples/nc_history_trace.exe).Path
Get-FileHash -Algorithm SHA256 $Trace
$T = 'interim/2026-11-22-impl-linear-solver'
$Out = 'work/nc96-compare-01'
python "$T/prepare_inputs.py" $Out --source 'D:/pc_runs/fine_20261122/nc96/input'
if ($LASTEXITCODE) { throw 'input preparation failed' }
```

`--source` は実際の `mesh.json`、`ground.json`、`stages.json` の所在に合わせる。過去の工程ファイルが `case.json` の場合は、初期入力であることを確認して新しい入力フォルダに `stages.json` としてバイトコピーする。スクリプトは元の3ファイルの SHA-256 を表示し、そのバイトを保存する。比較用の `stages-refined.json` は `settings.dense_lu_refinement=true` だけを変更する。

クラウドと同じ入力を作る場合は、新規出力名で `--source …` を `--cloud-original` に変更する。これは既存の `inputs('midpoint')` から Tet10 を一度だけ細分し、既往96要素と同じ自由度収容数750に設定したもの。既往の面集合別100 kPa手順の再生には `--cloud-face-sets` を使う。PC 原票との比較では同一の3入力を両方式に与える。

## 同一バイナリで既定と明示選択を比較

```powershell
$Dir = (Resolve-Path $Out).Path
foreach ($mode in @('baseline','refined')) {
  $caseFile = if ($mode -eq 'baseline') { 'stages.json' } else { 'stages-refined.json' }
  if (Test-Path (Join-Path $Dir $mode)) { throw "existing output $mode" }
  $clock = [Diagnostics.Stopwatch]::StartNew()
  $proc = Start-Process -FilePath $Trace -WorkingDirectory $Dir -PassThru -NoNewWindow `
    -ArgumentList @('direct','mesh.json','ground.json',$caseFile,$mode) `
    -RedirectStandardOutput "$Dir/$mode.stdout.log" -RedirectStandardError "$Dir/$mode.stderr.log"
  $peak = 0L
  while (-not $proc.WaitForExit(1000)) {
    $proc.Refresh()
    if ($proc.HasExited) { break }
    $peak = [Math]::Max($peak, $proc.PeakWorkingSet64)
    if ($clock.Elapsed.TotalHours -ge 6 -or $peak -ge 32000000000 -or $proc.PrivateMemorySize64 -ge 32000000000) {
      $proc.Kill()
      $proc.WaitForExit()
      throw "$mode stopped at resource budget; keep partial files"
    }
  }
  $proc.WaitForExit()
  "$mode seconds=$($clock.Elapsed.TotalSeconds) peak_working_set_bytes=$peak rc=$($proc.ExitCode)"
  if ($proc.ExitCode -ne 0) { throw "$mode process failed" }
}
```

各 `run.json` で `complete`、`failure_reason`、`last_committed.stage_progress`、最後の `attempts` の `progress_from/to` と `reason` を読む。受理数には geostatic.attempts の受理も含める。新方式の `linear_solves` では `iterations=1` が通常の LU 解、2 が1回の補正を経た解で、記録される `true_relative_residual` は元の行列積による。途中で棄却された Newton 試行の線形報告は Run に残らない場合がある。終了コード0だけでは全撤去も線形改善も判定できない。

## クラウドでの観測と保存行列の再測定

基準観測用の固定作業ツリーと実装用ツリーを分けた。以下は Linux シェルで、PR の実装用ツリーのルートから実行する。出力は新規の絶対パス、Rust 1.99.1 と Python 3.13、固定依存関係を使用する。観測器はソースの SHA をビルド前後で照合する。

```sh
set -eu
root=$(pwd)
out="$root/work/nc96-replay-01"
mkdir -p "$out"
base=518912cbf80b31e7a9b5736e338c8aa724681d77
git worktree add --detach "$out/fixed" "$base"
test "$(git -C "$out/fixed" rev-parse HEAD)" = "$base"
test -z "$(git -C "$out/fixed" status --porcelain)"
python "$out/fixed/tools/prepare_nc_linear_observer.py" --output "$out/observer"
T=interim/2026-11-22-impl-linear-solver
python "$T/prepare_inputs.py" "$out/inputs" --cloud-original
(
  cd "$out/inputs"
  ulimit -v 36000000
  NC_LINEAR_SYSTEMS="$out/rejected.jsonl" timeout 21600 "$out/observer/target/release/examples/nc_history_trace" direct mesh.json ground.json stages.json baseline
)
cargo build --locked --release --example nc_history_trace
(
  cd "$out/inputs"
  ulimit -v 36000000
  timeout 21600 "$root/target/release/examples/nc_history_trace" direct mesh.json ground.json stages-refined.json refined
)
```

保存済みの7件の96要素系と1件の12要素系は、非線形経路をもう一度解かずに検査できる。回帰試験は棄却と補正後受理を実装の元の行列積で検査する。Python の監査は保存した補正解と 80 桁積和、SVD、別の部分ピボット LU を使う。

```sh
F=tests/fixtures/linear-refinement
T=interim/2026-11-22-impl-linear-solver
cargo test --locked --release --lib stages::dense_refinement -- --nocapture
python "$T/audit_matrices.py" "$F/nc12.json" --solutions "$T/corrected-solutions.jsonl"
python "$T/audit_matrices.py" "$F/nc96.jsonl" --solutions "$T/corrected-solutions.jsonl"
```

12要素の原票は `interim/2026-11-08-nc-compensated-linear-replay/evidence.zip`（SHA-256 `4970434c69c08dfa30b8cd768d75ec80a8e7cf01b479e836a2903ffe99ed2bcc`）にある。`matrices.json` の `terminal-rejected:11` の列配列を行配列に転置し、`linear-results.json` の同一ラベルの基準解・残差を対応させた。96要素の `nc96.jsonl` は観測器の出力原文（SHA-256 `7e380e9067bcdfd1eaa77911705f1bec357d3e0a4cbdea782a8b338bc4047a90`）。別環境では演算順序に伴う経路差があり得る。進捗、目標、理由、行列 SHA を照合してから同一試行と判断する。

実装後の新しい終端も独立した観測ビルドで再実行した。補正した解の棄却行列は同じ監査器で再測定できる。

```sh
python "$T/audit_matrices.py" "$T/refined-terminal.jsonl"
```

Windows で P6-01 台帳の全原票を検証する場合は、CI と同様に最初のチェックアウトから Git の自動改行変換を無効にした作業コピーを使う。たとえば新規クローンには `git clone -c core.autocrlf=false --branch codex/plastic-spatial-convergence <PR のリポジトリ URL> <新規フォルダ>` を使う。今回の文書、保存行列、測定 JSON、追補台帳には個別の Git 属性も指定した。

## 資源見積もり

対象は96 Tet10 のみで、PC では1経路ずつ実行する。250節点の並進自由度の上限は750、750四方の binary64 行列は4.5 MB、行列と LU 因子の合計は約9 MB。撤去後の実測は348自由度。補正は同じ因子とベクトルを用いる。既往PCで538秒・49 MB、今回クラウドの既定経路で約538秒・最大観測32100 KiB、補正経路で約675秒・30380 KiB だった。PC では履歴・入出力も含め1経路1 GB・1時間程度を目安とし、経路が長引いた場合に備えて6時間・32 GBで打ち切る。停止した場合は部分結果と実測値を残す。768、6144要素は対象外。
