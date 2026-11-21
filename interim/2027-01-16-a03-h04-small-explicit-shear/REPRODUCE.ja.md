# 12 Tet10 の入力・封印・再実行

PR #13 のこの記録を含む checkout で、リポジトリルートを作業ディレクトリとする。Windows PC での再実行は任意。対象は高さ 0.4 m の界面、同じ 12 Tet10、54 節点である。子プロセスごとに 900 秒の上限を設け、各入力を順に解く。既存の12要素手順のメモリ目安はビルド 8 GB、求解・観測 12 GB 以内。Windows での実測はしていない。

必要なファイルは今回のフォルダ全体（全 part を含む）、01-07 明示せん断の Python と archive-index.json、01-06 不等厚の Python、12-30 A-03 追補と 11-26 A-03 の Python、src、examples、tools、Cargo.toml、Cargo.lock、rust-toolchain.toml、requirements*.txt と今回の回帰試験。Git の core.autocrlf=false で元のバイトを取得する。

## 保存値の照合

PowerShell から新しい出力先を使う。封印は全 part、結合 ZIP、全 member の長さと SHA-256 を照合して展開する。

~~~powershell
$ErrorActionPreference='Stop'
$D='interim/2027-01-16-a03-h04-small-explicit-shear'
$Old='interim/2027-01-07-a03-explicit-shear-initial-stress'
$Run=Join-Path 'work' ('a03h04-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Force $Run | Out-Null
$Py='.venv/Scripts/python.exe'
if (!(Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存の導入失敗' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
& $Py "$Old/a03shear_evidence.py" --verify $D --extract "$Run/sealed"
if ($LASTEXITCODE -ne 0) { throw '封印照合失敗' }
& $Py -m unittest discover -s tests -p test_zz_h04_small_shear.py -v
if ($LASTEXITCODE -ne 0) { throw '保存値の回帰失敗' }
& $Py "$D/a03h04_report.py" --run "$Run/sealed/formal" --output "$Run/from-seal"
if ($LASTEXITCODE -ne 0) { throw '再集計失敗' }
foreach ($file in @('measurements.json','signed-differences.jsonl')) {
    if ((Get-FileHash "$Run/from-seal/$file").Hash -ne (Get-FileHash "$Run/sealed/retained/report/$file").Hash) { throw "再集計差 $file" }
}
~~~

全受理場は各入力の *-fields.jsonl、順序の対応は *-layout.json、初期力は pre-solve-seed-forces.json、配置後の保持力は *-prepared.jsonl。measurements.json は共通点での全成分の生分母・絶対差・相対差と各経路の仕事・残差を含む。差が出た場合は出力と環境を保持する。

## 同じ入力の新規求解

~~~powershell
cargo build --locked --release --bin soil-fem-3d --example wet_birth_a03_probe
if ($LASTEXITCODE -ne 0) { throw '本体・既存観測 example のビルド失敗' }
& $Py "$D/a03h04_run.py" --output "$Run/formal"
if ($LASTEXITCODE -ne 0) { throw '求解未完了' }
& $Py "$D/a03h04_report.py" --run "$Run/formal" --output "$Run/reported"
if ($LASTEXITCODE -ne 0) { throw '観測未完了' }
~~~

完成した入力、途中状態、終了コード、報告を Linux 保存値と区別して残す。Windows 側の照合差を許容値で吸収しない。上限を超える一般空間確認は**規模上限のため確認不可・未達**。
