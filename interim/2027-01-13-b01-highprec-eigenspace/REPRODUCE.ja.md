# 保存行列の再実行

作業ディレクトリは PR #13 の本記録を含むリポジトリのルート。既存の証拠を読み、新規 ID の出力先に書く。対象は 25 節点・8 要素の 150 × 150 行列、支持後の照合はその主小行列 54 次。PC の操作実績はない。全体の受入試験や大規模な曲面解析を PC に引き渡す項目はない。

## 入力の照合と軽い回帰

Windows ではバイト照合用の checkout に最初から core.autocrlf=false を設定する。固定元コミットは 23c7b0ea48906723620b47eefb50e7481074d875。既存 checkout に一括の改行変換をかけない。入力 ZIP と目録は同コミットの tree から得た Git blob、ZIP 全体と各エントリーは元の目録の SHA-256 で検査する。欠落・不一致では停止する。

~~~powershell
$ErrorActionPreference = 'Stop'
$D = 'interim/2027-01-13-b01-highprec-eigenspace'
$Old = 'interim/2027-01-01-b01-curved-tri6-cloud-followup'
$Pin = '23c7b0ea48906723620b47eefb50e7481074d875'
git ls-tree $Pin -- "$Old/archive-index.json" "$Old/evidence.zip.part000"
if ($LASTEXITCODE -ne 0) { throw '固定 tree が必要' }
$Py = '.venv/Scripts/python.exe'
if (-not (Test-Path $Py)) {
  py -3.13 -m venv .venv
  if ($LASTEXITCODE -ne 0) { throw 'venv 作成失敗' }
}
& $Py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '固定依存の導入失敗' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
& $Py -m unittest discover -s tests -p 'test_zz_b01eig_spaces.py' -v
if ($LASTEXITCODE -ne 0) { throw '入力・封印・binary64 回帰の照合失敗' }
~~~

出力 ZIP の照合と展開は既存の読取専用検証器でも行える。封印には両精度の全固有値・全固有ベクトル、主角、残差、支持後の分解と実行ソースを収録した。各 JSON の行・列番号は節点順、各節点の並進 xyz、回転 xyz の順。

~~~powershell
$Run = Join-Path 'work' ('b01eig-' + [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffffffZ'))
if (Test-Path $Run) { throw '出力 ID が衝突' }
New-Item -ItemType Directory -Path $Run | Out-Null
& $Py 'interim/2026-12-31-a04-richards-wet-joint-cloud-followup/a04follow_evidence.py' --verify $D --extract "$Run/sealed"
if ($LASTEXITCODE -ne 0) { throw '封印照合失敗' }
~~~

## 任意の小型再計算

全 32 行列を 1 プロセスで順次実行できる。高精度演算には A-01 の先例と同じ mpmath 1.3.0 を明示導入する。通常の CI はこれを import しない。各コマンドは既存の出力先への再書込みを拒否する。失敗時はログと途中状態を保持する。

~~~powershell
& $Py -m pip install mpmath==1.3.0
if ($LASTEXITCODE -ne 0) { throw '任意精度依存の導入失敗' }
& $Py "$D/b01eig_solve.py" --workers 1 --output "$Run/raw"
if ($LASTEXITCODE -ne 0) { throw '固有値分解未完了' }
& $Py "$D/b01eig_measure.py" --raw "$Run/raw" --output "$Run/observations"
if ($LASTEXITCODE -ne 0) { throw '観測未完了' }
& $Py "$D/b01eig_supported.py" --output "$Run/supported.json"
if ($LASTEXITCODE -ne 0) { throw '支持部分空間の観測未完了' }
~~~

raw/ のベクトルは列。基底の符号・重複固有値内の基底選択は任意なので、再計算との比較には observations/ の主角・部分空間残差・射影量を用いる。SHA 固定は封印した生値そのものを識別する。binary64 の LAPACK が別実装の環境で任意基底のバイト一致を要求する検査は行わない。

Linux でも同じ引数で実行できる。最大回転差が生じた 1 成分の符号付き内訳は supported.json の ordinary_delta_rad、saved_decimal80_delta_rad、low9_error_delta_rad、other45_error_delta_rad にある。各形状の identity 以外の出力は、identity 基底を保存回転行列で移した後の主角・全固有値差を含む。

## 停止記録 ZIP の初回固定

interruption-record.zip には以前の封印値がなく、今回が初回固定。2026-11-19 に作業ツリーと Git インデックスのバイト一致および ZIP CRC を検査し、SHA-256 `727317c3cc3171b0423e489e76123cd48b377321626134a5b7e92d88aabb40ff` を archive-index.json の supplemental_archives に追記した。
