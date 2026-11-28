# NC-01b の再実行（Windows PC）

リポジトリのルートから実行する。対象は PR #13 の NC-01b 登録を含むコミット。Python はリポジトリの 3.13 環境と `requirements-reference.lock.txt`、Rust は `rust-toolchain.toml` と `Cargo.lock` の指定を使う。ここに記録した測定はクラウド Linux のもの。Windows PC での計測は未実施。

## 規模と打切り

|要素数|節点・並進成分上限|本体のクラウド全経路実測|PC の計画値（推定）|
|---|---|---|---|
|12 Tet10|45 節点、135 成分|0.85～1.26 秒|1 経路 10 分・4 GiB を執行上限とする。|
|96 Tet10|225 節点、675 成分|10.54～12.97 秒|1 経路 30 分・4 GiB を執行上限とする。|

96 の全成分を密に置いた倍精度行列 1 枚は 3,645,000 bytes。実際は元の疎 GMRES-ILU0 を用いる。行列以外に材料履歴、フレーム、プロセス本体が必要で、全 RSS はこの行列サイズからは求まらない。Linux の本計算はアドレス空間 4 GiB の制限内で完了した（ピーク RSS 未測定）。Windows の既存ラッパーは子を停止状態で Job Object に入れ、プロセス群のコミットメモリと時間を制限し、ピーク値と終了理由を保存する。全経路を逐次実行し、各上限を超えた場合は未達として終了する。途中再開や自動の試行追加は行わない。生成器の選択肢は 12 と 96。768・6144 要素の実行は対象外。

## 入力生成と全撤去

以下は 96 要素の全3経路。12 要素には `$Elements = 12` を使う。出力先は新規フォルダを指定する。旧 PC 転記の `nc96` フォルダと区別する。

```powershell
$ErrorActionPreference = 'Stop'
git merge-base --is-ancestor ee68ee2ef8e904fd71643f8ed566d2e49eeaebd0 HEAD
if ($LASTEXITCODE -ne 0) { throw 'NC-01b の入力を含む版が必要' }
$Revision = git rev-parse HEAD
$Py = (Resolve-Path .venv/Scripts/python.exe).Path
$env:PYTHONUTF8 = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'release build failed' }
$Solver = (Resolve-Path target/release/soil-fem-3d.exe).Path
$Guard = 'interim/2027-01-08-pc-failure-diagnosis/pcfail_run.py'
$Study = 'interim/2026-11-22-impl-nc01-option-b-physical-input'
$Elements = 96
$Seconds = if ($Elements -eq 12) { 600 } else { 1800 }
$Run = "runs/NEW-nc01b-$Elements"
if (Test-Path $Run) { throw '未使用の出力先を指定する' }
New-Item -ItemType Directory $Run | Out-Null
$Revision | Set-Content -Encoding ascii "$Run/revision.txt"
foreach ($Profile in @('k0','roller','pressure')) {
    $Dir = "$Run/$Profile"
    & $Py tools/nc01b_inputs.py --profile $Profile --elements $Elements --output $Dir
    if ($LASTEXITCODE -ne 0) { throw "入力または施工前の釣合いを確認 $Profile" }
    & $Py $Guard --output "$Dir/process" --timeout $Seconds --memory-gib 4 -- $Solver solve-stages "$Dir/mesh.json" "$Dir/ground.json" "$Dir/stages.json" "$Dir/native.json"
    $SolveCode = $LASTEXITCODE
    Get-Content "$Dir/process/process.json"
    if (Test-Path "$Dir/native.json") {
        & $Py "$Study/summarize.py" "$Dir/native.json" "$Dir/initial-balance.json"
        if ($LASTEXITCODE -ne 0) { throw "受理履歴の照合失敗 $Profile" }
    }
    if ($SolveCode -ne 0) { throw "経路停止 $Profile。出力を保持し process.log と failure_reason を確認" }
}
```

`initial-balance.json` は反力・合力を含む別積分。入力3ファイルと `native.json` の SHA-256、`revision.txt`、`process.json`、`process.log` を同じ経路に保管する。終了コード、初期釣合い、除荷完了時の力残差、撤去開始のゼロ解放残差、受理した撤去 λ、最終の活動要素数・保持荷重、最大の真の線形残差を読む。時間打切りで本体結果が保存されなければ、最後の受理 λ は不明と記録する。チェックポイントの `stage_progress` は初期化・載荷・除荷中にも用いるため、`excavate_upper` の受理増分と照合する。

## 保存済み結果の再監査

```powershell
$Evidence = "runs/NEW-nc01b-cloud-evidence"
if (Test-Path $Evidence) { throw '新規展開先が必要' }
(Get-FileHash -Algorithm SHA256 "$Study/cloud-runs.zip").Hash
Expand-Archive "$Study/cloud-runs.zip" $Evidence -ErrorAction Stop
& $Py "$Study/summarize.py" "$Evidence/96-pressure/run.json" "$Evidence/96-pressure/initial-balance.json"
& $Py interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
if ($LASTEXITCODE -ne 0) { throw 'P6-01 証拠照合失敗' }
```

ZIP の期待 SHA と P6-01 の固定コミット・blob・台帳 SHA は [README.ja.md](README.ja.md) に掲載。ZIP 内の経路別ファイル SHA は [cloud-runs.json](cloud-runs.json) にある。本体ファイルは実行時間・バイナリ識別子を含むため、新しい PC 実行では入力と物理量・残差・到達状態を比較する。

元の NC-01 は未達・分類 D。
