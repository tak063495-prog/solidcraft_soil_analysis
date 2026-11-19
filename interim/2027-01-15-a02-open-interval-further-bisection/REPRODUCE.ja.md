# 同じ 24 Tet10 の結果を読む・再計算する

PR #13 の本記録を含む checkout、Python 3.13 と Rust 1.99.1 を使用する。Git は core.autocrlf=false で checkout する。sparse の場合は新フォルダ全体（全 evidence.zip.part*）、01-05 A-02 の *.py、12-29・12-21・11-25 の A-02 スクリプトと 12-29 の a02follow_service.rs、Cargo.toml、Cargo.lock、rust-toolchain.toml、src、tools、examples、requirements*.txt、今回の tests/test_zz_further_a02_bisection.py が必要。01-05 の対照 L0 は本封印に元 member と同一のバイトで収録した。

## 封印の照合と表の再作成

リポジトリルートの PowerShell で実行する。出力先は実行ごとに新規作成する。

~~~powershell
$ErrorActionPreference = 'Stop'
$D = 'interim/2027-01-15-a02-open-interval-further-bisection'
$Old = 'interim/2027-01-05-a02-n8-contact-cap-open-interval'
$Run = Join-Path 'work' ('a02bisect-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Force -Path $Run | Out-Null
$Py = '.venv/Scripts/python.exe'
if (!(Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存の導入失敗' }
$env:CARGO_BUILD_JOBS='2'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
& $Py "$Old/a02n8_evidence.py" --verify $D --extract "$Run/sealed"
if ($LASTEXITCODE -ne 0) { throw '封印照合失敗' }
& $Py -m unittest discover -s tests -p test_zz_further_a02_bisection.py -v
if ($LASTEXITCODE -ne 0) { throw '保存値の回帰失敗' }
& $Py "$D/a02bisect_report.py" --paths "$Run/sealed/formal" --output "$Run/reported"
if ($LASTEXITCODE -ne 0) { throw '観測再集計失敗' }
foreach ($name in @('first-events.csv','interval-differences.csv','measurements.json')) {
    if ((Get-FileHash "$Run/reported/$name").Hash -ne (Get-FileHash "$D/$name").Hash) { throw "保存値と異なる $name" }
}
~~~

## 同じ入力での任意の Windows 求解

対象は 81 節点、24 Tet10、32 接触 GP、通常接触上限 0.05 m、M=4、N=8 のみ。3 水準それぞれ本体・独立を順次解く。本体は観測なしの whole も各回に起動し、計 9 経路になる。01-05 のビルド器は指定の別フォルダにソースをコピーし、そのコピーから専用 probe を作成する。旧手順の小型ケース予算（ビルド 8 GB 以内・主経路 16 分以内）を参照した。Windows での時間・ピークメモリは実測していない。新規求解では各 OS の生値を別保存する。

~~~powershell
& $Py "$Old/a02n8_build.py" --output "$Run/build"
if ($LASTEXITCODE -ne 0) { throw '隔離ビルド失敗' }
New-Item -ItemType Directory -Path "$Run/formal" | Out-Null
Copy-Item -Recurse "$Run/sealed/formal/L0" "$Run/formal/L0"
foreach ($level in 1..3) {
    New-Item -ItemType Directory -Path "$Run/formal/L$level" | Out-Null
    foreach ($engine in @('native','independent')) {
        & $Py "$D/a02bisect_run.py" --level $level --engine $engine --build "$Run/build" --output "$Run/formal/L$level/$engine"
        if ($LASTEXITCODE -ne 0) { throw "経路未完了 L$level/$engine" }
    }
}
& $Py "$D/a02bisect_report.py" --paths "$Run/formal" --output "$Run/new-report"
if ($LASTEXITCODE -ne 0) { throw '新規求解の観測未完了' }
~~~

各段の全 GP の初回受理は first-events.csv、各対応区間の生差ノルムは interval-differences.csv、符号付きの全成分は signed-differences.json、元の終点場は各経路の *-fields.jsonl にある。GPD WIN Max 2（利用可能メモリ約 48 GB）での任意確認はこの小型再現まで。規模上限を超える一般確認は**規模上限のため確認不可・未達**。
