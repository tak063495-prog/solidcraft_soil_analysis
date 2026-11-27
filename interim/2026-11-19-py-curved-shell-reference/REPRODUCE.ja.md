# 再現手順（Windows PC中心）

保存先ラベル2026-11-19、参照系列・新規試験2026-11-11、既存回帰試験2026-11-12（日本時間）。B-01は未達。大きいメッシュ、実CAD、GUI/GPUはクラウド未実行。本参照のためにRust・既存Pythonソルバは不要。既存Python回帰試験だけは従来どおりCLIを使う。

## 準備

Windows PowerShell、Python3.13、64bitを想定。リポジトリのルートで実行する。新規出力先を毎回変え、保存証拠を上書きしない。

```powershell
git fetch origin codex/py-curved-shell-reference
git switch codex/py-curved-shell-reference
py -3.13 -m venv work/curved-shell-venv
$py = (Resolve-Path work/curved-shell-venv/Scripts/python.exe).Path
& $py -m pip install -r tools/py_curved_shell_reference/requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
& $py tools/py_curved_shell_reference/check_reference.py
if ($LASTEXITCODE -ne 0) { throw 'standalone checks failed' }
& $py tools/py_curved_shell_reference/reference.py --levels 4 8 16 32 --out work/curved-shell-pc-small
if ($LASTEXITCODE -ne 0) { throw 'reference run failed' }
```

27試験・skip0が必要。小系列の目安は普通のPCで10–60秒（今回Linux実測は求解合計約6秒）。出力は `results.json` と3例×4水準の全節点CSV、標準出力に各水準の実測。`results.json` が作られない途中終了を成功にしない。依存インストールの時間は別。

Linuxで同じ新規検証を行う場合：

```bash
python3 -m venv work/curved-shell-venv
work/curved-shell-venv/bin/python -m pip install -r tools/py_curved_shell_reference/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 work/curved-shell-venv/bin/python tools/py_curved_shell_reference/check_reference.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 work/curved-shell-venv/bin/python tools/py_curved_shell_reference/reference.py --levels 4 8 16 32 --out work/curved-shell-linux-small
```

## PC専用の大きいh系列

クラウドでは最大n=32まで実測。次の64/128は**未実行**で、Windows PC（推奨32GB RAM、SSD）だけで実行する。時間は移植後の予算目安であり測定値ではない。

```powershell
& $py tools/py_curved_shell_reference/reference.py --levels 16 32 64 --allow-large --out work/curved-shell-pc-64
if ($LASTEXITCODE -ne 0) { throw 'n64 series incomplete or failed' }
# n64の経過時間とRAMを確認してから別フォルダーへ進む
& $py tools/py_curved_shell_reference/reference.py --levels 32 64 128 --allow-large --out work/curved-shell-pc-128
if ($LASTEXITCODE -ne 0) { throw 'n128 series incomplete or failed' }
```

n64は4,225節点/21,125自由度、n128は16,641節点/83,205自由度。各コマンドは全3例を求解する。目安はn64系列1–10分、n128系列5–30分、ピークRAM数GB以上。SciPyの疎直接分解・Python組立リストでメモリーが増える。実測で余裕がなければ中止し、完了/不完了、経過時間、例名、n、例外を保存する。終了・OOM・残差上限超過・中断を合格にしない。クラウドでn>32を試さない。`--allow-large`なしではn>32を拒否する。

**判定**：内部検査は、有限な解、自由自由度の相対方程式残差≤1e−7、n=2/4の無拘束零モード6・余剰0・有意負固有値なし・解析剛体残差<1e−9、拘束後固有値>1e−9。これらは線形求解・数値ランクの検査であり、文献変位の受入基準ではない。文献誤差と観測次数は全部保存し、負の次数や再増加を消さない。変位許容値・最低次数は**未定義**なので、B-01合格の判定はできない。

MRの点荷重直下変位には特異性の留保がある。大きいnの円筒・半球が従来文献値を越えて増大しても、荷重分散や参照値の変更で合わせない。正則化を別研究で行う場合は別問題・別データとして扱い、この入力・値を変更しない。

## 将来の本体曲面Tri6出力との比較

CSVは `node_id,x,y,z,ux,uy,uz`（node_idは任意）、JSONは同じキーのオブジェクト配列、または次の形式。座標・並進はベンチマークJSONと同じ単位・グローバル系。IDの連番や本参照と同じメッシュは不要。ただし評価物理点の節点がちょうど1個必要。

```json
{"benchmark_id":"hemisphere_hole","nodes":[{"node_id":23,"x":10,"y":0,"z":0,"ux":0.093,"uy":0,"uz":0}]}
```

本体は現時点で曲面Tri6を拒否するため、本体比較用の実出力はまだない。次は本参照自身のCSVを使った形式確認であり、本体との比較実績ではない。

```powershell
& $py tools/py_curved_shell_reference/compare.py --benchmark scordelis_lo --input work/curved-shell-pc-small/scordelis_lo-n32.csv --out work/compare-roof.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined tolerance, exit 2' }
& $py tools/py_curved_shell_reference/compare.py --benchmark pinched_cylinder --input work/curved-shell-pc-small/pinched_cylinder-n32.csv --out work/compare-cylinder.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined tolerance, exit 2' }
& $py tools/py_curved_shell_reference/compare.py --benchmark hemisphere_hole --input work/curved-shell-pc-small/hemisphere_hole-n32.csv --out work/compare-hemisphere.json
if ($LASTEXITCODE -ne 2) { throw 'expected undefined tolerance, exit 2' }
# 将来は --input をユーザーPC上の本体CSV/JSONパスに置き換える
```

出力の差・status・passed・許容値根拠を見る。現状は正常な入力でも `status="未定義" / passed=null / relative_tolerance=null`。終了2を成功0へ変換してCI合格にしない。評価点欠損/重複、異なるbenchmark_id、NaN/Inf等は入力不正で終了3。比較1成分だけで全場やB-01を受け入れない。

## 既存Python全件（設定変更なし）

既存回帰試験は別の従来環境で実行する。独立参照自体はこのCLIを一切呼ばない。従来依存のセットアップ後、Rust1.99.1でrelease CLIをビルドし、既存4群すべてを実行する。変更対象は追加ファイルだけなのでCI設定・収集件数条件は編集しない。

```powershell
py -3.13 -m venv work/regression-venv
$regpy = (Resolve-Path work/regression-venv/Scripts/python.exe).Path
& $regpy -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'regression dependencies failed' }
rustup toolchain install 1.99.1 --profile minimal
if ($LASTEXITCODE -ne 0) { throw 'Rust toolchain setup failed' }
cargo +1.99.1 build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'release build failed' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
    & $regpy tools/run_ci_python_tests.py --group $group --output "work/curved-shell-regression/$group"
    if ($LASTEXITCODE -ne 0) { throw "regression failed: $group" }
}
```

Linuxも `python tools/run_ci_python_tests.py --group <群> --output <新規フォルダー>` を4群実行する。実行全件・失敗0・エラー0・skip0・expected failure0・unexpected success0と、4群のIDの和集合が既存収集全体に重複なく一致することが必要。最終の既存件数は532（265＋259＋4＋4）。保存証拠はGB規模で、全件回帰には殻参照とは別に数十分程度かかる（初回取得時間は回線依存で別）。今回の平均応力証拠群だけで約793秒だった。件数・ログは `validation.json` と照合する。新規27試験は上の明示コマンドで別に実行する。

Windows CIのhead SHAと各群の実行完了を別途確認する。タイムアウト・キャンセル・旧SHAの成功を新先端の成功へ流用しない。数値上の未達は回帰成功と独立に保持する。
