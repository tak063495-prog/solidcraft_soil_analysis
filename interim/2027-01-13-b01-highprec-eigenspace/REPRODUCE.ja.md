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

## Linux での NC 履歴照合の切り分け

比較には `f460f562fa42ebfe72f11cd1e3cdbbb07a241762` と `f74e09a5288f2d126894d40c322fb3e446946aae` の別々の checkout を使った。それぞれで NC の保存アーカイブ 328 part を用意し、同じ固定依存の Python と単一スレッドの BLAS を使う。以下の `BASE` と `HEAD_TREE` はこの順の checkout、`OBSERVER` は本フォルダの `b01eig_nc_difference.py` の絶対パス、`OUT` は新規の絶対出力先とする。二つの群は単独で順に実行し、それぞれの終了コードと標準エラーを保存する。

~~~bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
(cd "$BASE" && python tools/run_ci_python_evidence_shards.py --group python-mean-evidence --output "$OUT/base-mean") > "$OUT/base-mean.log" 2>&1
base_exit=$?
(cd "$HEAD_TREE" && python tools/run_ci_python_evidence_shards.py --group python-mean-evidence --output "$OUT/head-mean") > "$OUT/head-mean.log" 2>&1
head_exit=$?
printf 'base=%s head=%s\n' "$base_exit" "$head_exit"
~~~

実際の単独実行は旧先端で `Ran 0 tests in 721.750s`、現先端で `Ran 0 tests in 715.360s`。両方で `FAILED (errors=1)`、`ValueError: saved mean-stress/complement partition differs from own history` を出した。両コミットで NC ツリー・試験・依存ロックの Git オブジェクトは一致する。[ログ・Git オブジェクト比較](linux-mean-comparison-interrupted.zip)

成分差は現先端の共通封印から 1 回だけ測定した。前回停止したスクリプトでは `sources/examples/` の入力が欠けていたため、元の目録からメッシュ・材料・ステージを読み、各バイト数と SHA-256 を照合してから展開する。既存の失敗は上記 ZIP 内の base-observe.log に残した。

~~~bash
python -I "$OBSERVER" --root "$HEAD_TREE" --output "$OUT/head-partition.json"
~~~

34 ケースのうち 29 ケース、483 成分に差があり、最大絶対差は 7.105427357601002e−15 kN（節点 19、z 方向、補成分）。出力の `node_zero_based` と `dof` は 0 始まり。符号付き差は再計算値から保存値を引いた値で、出力に両値の binary64 十六進表現を含む。試験結果は「f460f56 でも再現する Linux 環境依存の既存失敗（B-01 と無関係、原因未特定）」として扱う。[観測した全差分](nc-partition-differences.json)

## spatial 群の旧先端対照と確定件数

旧先端 f460f56 の spatial 群を同じ Linux 環境で 1 回実行した。出力は `Ran 0 tests in 169.207s`。現先端 f74e09a5 の保存ログは `Ran 0 tests in 179.791s`。両群は `FAILED (errors=1)`、終了 1。旧先端の群では採取用ハードリンクが `OSError: [Errno 18] Invalid cross-device link` となり、照合行に到達しなかった。JSON 差分採取時に旧先端の未変更の setUpClass 内の呼出し（remeasure.run、Linux の元引数）を保存先指定で直接実行すると、現先端の群と同じ `ValueError: remeasured JSON differs` が同じ行で再現した。群自体は再実行していない。NC 空間証拠ツリーは両方とも `2f093479f0c46e194e14f007e2616f51c276d2cd`。mean と spatial はいずれも Linux 環境依存の既存失敗（B-01 と無関係、原因未特定）として記録した。

旧先端で直接実行した元の照合が生成した JSON を保持し、保存 JSON と 1 回照合した。非一致の数値キーは 49 件。絶対差が最大のキーは `/paths/native-left256/remaining_holding_norm_kn` で、保存値 `3.883669499505146`、再測定値 `3.8836694995051464`、再測定 − 保存 `+4.4408920985006262e-16`。これは JSON に記録された各値の単位での数値差。異なる物理量を横断した誤差判定には使用しない。

**未実行 8 件・setUpClass エラー 2 件（既存）**。Rust は 467/467 成功。Linux Python は 998 件中 988 成功、PowerShell の skip 2。旧先端での対照実行は本数に重ねて加えていない。B-01 は未達、`passed=false`。[両コミットの出力、再測定 JSON、全差分キーと実行スクリプト](linux-spatial-comparison.zip) ／ [確定した Linux 検証記録](validation.json)

旧先端の群の元出力は base-spatial.log。capture-spatial.py の一時 JSON 採取処理で発生したハードリンクエラーも保存した。照合そのものを再現するには、旧先端 checkout のルートで固定依存・単一 BLAS スレッドを使用し、同梱の direct-spatial.py に新規出力先の絶対パスを渡す。元の Linux 引数で remeasure.run を呼び、終了 1 の後もそこに measurements.json が残る。compare-spatial.py に旧先端の interim/2026-11-10-nc-spatial-release/measurements.json、出力先の measurements.json、新規の差分出力パスを順に渡す。base-spatial-direct.log に元の照合の例外全文、spatial-json-differences.json に両 SHA-256 と値の異なる JSON Pointer 全件を収録した。
