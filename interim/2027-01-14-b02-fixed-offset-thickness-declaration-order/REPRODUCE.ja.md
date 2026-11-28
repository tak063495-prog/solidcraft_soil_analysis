# 12/96 Tet10 の再実行

PR #13 の本記録を含む checkout を使用する。Windows は GPD WIN Max 2、使用可能メモリ約 48 GB を想定。以下は利用者が任意に再現する手順。作業 ID ごとに新しい出力先を作る。実行範囲は 12/96 Tet10、最大 225 土節点、8 枚の Tri6 と 1 本の梁。各子プロセスの既存の制限は 7200 秒で、1 日を超える経路は対象外。従来の同型ケースのメモリ予算は求解 1〜4 GB、ビルド 2〜8 GB。Windows での今回の全経路の時間・ピークメモリは未計測。

## 保存値の照合

Git の改行変換を無効にした checkout と、Python 3.13、Rust 1.99.1 を使う。sparse checkout の場合は新フォルダ全体と旧 B-02 のヘルパー、tools、examples、src を用意する。分割 ZIP は全 part が必要。

~~~powershell
$ErrorActionPreference = 'Stop'
$D = 'interim/2027-01-14-b02-fixed-offset-thickness-declaration-order'
$Previous = 'interim/2027-01-02-b02-offset-composite-cloud-followup'
$Old = 'interim/2026-12-01-tri6-offset-composite-b02'
$Run = Join-Path 'work' ('b02split-' + [guid]::NewGuid().ToString('N'))
New-Item -Type Directory -Path $Run -Force | Out-Null
$Py = '.venv/Scripts/python.exe'
if (!(Test-Path $Py)) { py -3.13 -m venv .venv }
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存を導入できませんでした' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
& $Py "$Previous/b02follow_evidence.py" --verify $D
if ($LASTEXITCODE -ne 0) { throw '封印照合失敗' }
& $Py -m unittest discover -s tests -p test_zz_b02split_controls.py -v
if ($LASTEXITCODE -ne 0) { throw '小型回帰失敗' }
~~~

封印の展開には検証コマンドに `--extract "$Run/sealed"` を加える。各ケースの mesh・ground・stages、受理状態列、比較、状態ごとの仕事と残差、コマンド・終了コードを読める。元ゲートの相対残差は `original_relative`。同じ成分を固定順で再集計した値は `relative_fixed_order`。

## 同じ対照を新しく求解

~~~powershell
& $Py "$Old/build_observer.py" "$Run/build"
if ($LASTEXITCODE -ne 0) { throw '観測器ビルド失敗' }
& $Py "$D/b02split_run.py" --observer target/release/b02-observer.exe --output "$Run/controls"
if ($LASTEXITCODE -ne 0) { throw '経路未完了。ログと部分出力を保持' }
& $Py "$D/b02split_report.py" --run "$Run/controls" --output "$Run/measurements.json"
if ($LASTEXITCODE -ne 0) { throw '観測集計失敗' }
~~~

Linux では同じ引数で仮想環境の Python と拡張子なしの観測器を指定する。各ケース名の t は m。12/96 それぞれに既定順 3 板厚、薄板で要素だけ・面だけ・板だけの反転 3 条件がある。腕長と物理節点の番号は全ケースで同一。各実装の積分、荷重履歴、支持、施工順を保持する。再実行結果は元の封印に上書きせず、実行環境を添えて比較する。

PC で任意に確認するのはこの小型ケースの OS 間再現まで。768・6144 要素やメモリ 40 GB 超・1 経路 1 日超が必要な一般確認は**規模上限のため確認不可・未達**。
