# 乾燥 JOINT 第4反復の再現

PR #13 の本記録を含む checkout を使う。GPD WIN Max 2（利用可能メモリ約48 GB）での任意手順。24 Tet10・32接触点・153自由度の密行列2経路を順次計算する。ビルドは2ジョブ、BLAS は1スレッド。Windows PC の占有メモリ・所要時間は実測していない。一般空間収束は規模上限のため確認不可・未達。

## 保存値と小さい再計算

リポジトリ直下の PowerShell で実行する。初回 checkout 前に core.autocrlf=false を設定する。既存の固定 requirements、Rust toolchain、tools、examples、src、A-02 の 11-25、A-04 の 12-31 と 01-12 の参照スクリプト、このフォルダと新規試験が必要。今回の封印は1 part、29 member、4,521,093 bytes。結合 ZIP の SHA-256 は `71a7aa6ee8686dc07846744202a34c7095bed9696b6179adc5eae0a3b54880a9`。

~~~powershell
$ErrorActionPreference='Stop'
$D='interim/2027-01-17-a04-iteration4-tangent-rhs-correction'
$Verify='interim/2026-12-31-a04-richards-wet-joint-cloud-followup/a04follow_evidence.py'
$Run=Join-Path 'work' ('a04it4-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory $Run | Out-Null
$Py='.venv/Scripts/python.exe'
if (!(Test-Path $Py)) {
    py -3.13 -m venv .venv
    if ($LASTEXITCODE -ne 0) { throw 'venv 作成失敗' }
}
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の導入失敗' }
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'
& $Py $Verify --verify $D --extract "$Run/seal"
if ($LASTEXITCODE -ne 0) { throw '封印の長さまたは SHA が不一致' }
& $Py -m unittest discover -s tests -p test_zzz_a04it4_tangent.py -v
if ($LASTEXITCODE -ne 0) { throw '保存成分の軽量照合失敗' }
~~~

`formal/{000,010}/native.jsonl` が本体の実履歴。linear 行は全接線、右辺、線形解、物理自由度の対応を含む。accepted 行は第4補正の受理。independent-path.json は独立参照の実履歴（0〜2）、independent-fixed.json は本体の第4状態で再評価した比較用の系。両側とも右辺は外力−内力、変位は m、節点力は kN、接線は kN/m。retained/a04it4-final40 と final55 の数値配列は十進文字列。各補正の順序は ids に記載した物理自由度順である。

## 保存系から再集計

~~~powershell
foreach ($Digits in @(40,55)) {
    & $Py "$D/a04it4_report.py" --raw "$Run/seal/formal" --digits $Digits --output "$Run/recomputed-$Digits"
    if ($LASTEXITCODE -ne 0) { throw "十進 $Digits 桁の再計算失敗" }
    foreach ($File in @('000.json','010.json','elastic.json','orders.json','summary.json')) {
        $Saved=Get-FileHash "$Run/seal/retained/a04it4-final$Digits/$File"
        $Fresh=Get-FileHash "$Run/recomputed-$Digits/$File"
        if ($Fresh.Hash -ne $Saved.Hash) { throw "再計算の差 $Digits/$File" }
    }
}
~~~

## 入力からの再実行

~~~powershell
& $Py "$D/a04it4_build.py" --output "$Run/build"
if ($LASTEXITCODE -ne 0) { throw '隔離プローブのビルド失敗' }
& $Py "$D/a04it4_run.py" --probe "$Run/build/target/release/a04it4_dry.exe" --output "$Run/fresh"
if ($LASTEXITCODE -ne 0) { throw '乾燥対照または同一状態照合の失敗' }
& $Py "$D/a04it4_report.py" --raw "$Run/fresh" --output "$Run/fresh-report"
if ($LASTEXITCODE -ne 0) { throw '新規求解の比較失敗' }
~~~

同じ入力 SHA、順・接触状態・停止反復・生の成分差を新しい ID に保存する。環境が異なる新規求解は、保存した系を同じ順で再集計する手順と区別する。失敗時は出力・エラー・途中状態をそのまま残す。旧 PC の実行や失敗記録を再実行先にしない。
