# Windows PCでの独立参照・比較の再現

起点は `2d2d54c48070581d088959b52984ceb348e0e85e`、ブランチは
`codex/py-transient-water-reference`、PRのbaseは `codex/plastic-spatial-convergence`。
新規ファイルだけ。C-04/A-04・NC-01・P6は未達を維持する。
以下はリポジトリルートでの PowerShell 7 用コマンド。
既存ソース・テスト・CI・既存interimを編集せず、新しい `work/` の出力へ記録する。
ディレクトリ名にある日付は保存先ラベルであり、実施日ではない。

## 最小環境と小規模確認

CPython 3.13、NumPy 2.5.4、SciPy 1.18.1のみ。コンパイラや本体CLIは不要。
依存導入は回線次第で1〜5分、検査・smoke系列は通常各5〜30秒を見込む
（Linux実測は各約3.5秒。Windows時間は未測定）。

```powershell
$ref = 'interim/2026-11-19-py-transient-water-reference'
$out = 'work/py-water-reference-pc-01'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
New-Item -ItemType Directory -Force $out | Out-Null
py -3.13 -m venv work/py-water-venv
if ($LASTEXITCODE -ne 0) { throw 'venv failed' }
$python = (Resolve-Path 'work/py-water-venv/Scripts/python.exe').Path
& $python -m pip install -r "$ref/requirements.lock.txt"
if ($LASTEXITCODE -ne 0) { throw 'install failed' }
& $python -m pip freeze | Out-File -Encoding utf8 "$out/packages.txt"
& $python "$ref/checks.py" -v 2>&1 | Tee-Object "$out/new-checks.log"
if ($LASTEXITCODE -ne 0) { throw 'independent checks failed' }
& $python "$ref/diagnose.py" --profile smoke --output "$out/smoke.json"
if ($LASTEXITCODE -ne 0) { throw 'reference refinement checks failed' }
```

合否：新規24件が全て成功し、失敗・エラー・skipが0。
`smoke.json` の `reference_checks_passed=true` と、各ケースの時間次数・空間次数・
時間制御差・質量収支を確認する。各検査の数値条件はREADMEの「刻みの分離」に記載。
計算が中断されJSONがない/途中、終了コードが非ゼロなら、成功に数えない。
参照実装自身が通っても本体や湿潤JOINTの受入にはならない。

## PC専用の大きめの系列

クラウドでは未実行。時間方向は160セル固定、40/80/160/320刻み、
空間方向は20/40/80/160セルで高精度時間積分とその100倍厳しい対照を使う。
通常30秒〜3分程度の見込み（Windows未測定、PC性能に依存）。
これでも一般3Dや長時間Richardsの実務計算とは異なる小型1D参照である。

```powershell
& $python "$ref/diagnose.py" --profile windows-local --output "$out/windows-local.json" 2>&1 |
  Tee-Object "$out/windows-local.log"
if ($LASTEXITCODE -ne 0) { throw 'local refinement did not pass' }
```

合否条件はsmokeと同じ参照自己検査。観測次数・時間制御差・収支を全水準で保存する。
有限の最細格子を厳密解と呼ばず、解析式との差と半離散の時間誤差を区別する。
結果が想定条件を満たさないときは元のJSONを保存し、許容値や参照式を合わせ直さない。

## CSV/JSON比較器の人工例

次は解析式から作るアダプターの例で、本体出力でも独立数値求解の証拠でもない。
JSONとCSV＋sidecarの両方を生成する。約1秒。

```powershell
& $python "$ref/compare.py" --demo "$out/demo"
if ($LASTEXITCODE -ne 0) { throw 'demo generation failed' }
& $python "$ref/compare.py" "$out/demo/terzaghi.json" --output "$out/demo-json-comparison.json"
if ($LASTEXITCODE -ne 2) { throw 'expected valid diagnostic with undefined acceptance (exit 2)' }
& $python "$ref/compare.py" "$out/demo/richards.csv" `
  --metadata "$out/demo/richards.metadata.json" --output "$out/demo-csv-comparison.json"
if ($LASTEXITCODE -ne 2) { throw 'expected diagnostic exit 2' }
```

比較器の終了コード：`--demo` 生成0、正常な計測だが**受入未定義2**、不正入力1。
引数不備もargparseによって2になるため、診断JSONの存在・schema・`accepted=null`を
必ず併せて確認する。終了コード2だけで「正常計測」や「合格」と判断しない。
数値誤差がゼロの人工例も `accepted=null`、C-04/A-04は未達である。

## 将来の本体出力からの入力契約

この追加は本体のCSV/JSON出力形式を変更しない。元の出力を保持し、PC上の別アダプターで
**同じベンチマークの確定時刻・物理点**を下記へ変換する。現時点で本体に指数型Richards
モデルがないなら、VG/Mualemを `richards` とラベル付けして使わない。
一般3D・JOINTの場をこの1D参照へ勝手に縮約しない。

CSVは `time,x,head,flux` の4列。列値は有限数、時間は少なくとも2時刻、
同一(time,x)の重複は禁止。JSONは同じ4キーの辞書を `samples` 配列に置く。
座標は0〜1 m、比較時刻はTerzaghi 0.02〜0.2 h、変動水位0〜0.2 h、Richards 0.5〜5 h。
境界を含め任意の物理点・時刻で直接解析式を評価し、自動補間しない。
比較窓外の初期特異点を無断で除いた時系列の合格ではなく、明示した窓の診断である。

メタデータはJSONでは最上位、CSVでは `--metadata` のJSONに置く。必須値は
`reference.metadata(case)` と一致させる。

```json
{
  "schema": "transient-water-1",
  "case": "richards",
  "benchmark": "fixed-v1",
  "length_unit": "m",
  "time_unit": "h",
  "head_unit": "m",
  "flux_unit": "m/h",
  "flux_kind": "instantaneous",
  "head_kind": "total_head",
  "coordinate": "up_from_water_table",
  "flux_sign": "positive_increasing_coordinate",
  "complete": true,
  "run_status": "completed",
  "samples": []
}
```

上の空配列は**形式説明のみで単独では拒否される**。有効な人工時系列例は `--demo` で生成する。
`provenance` 等に元出力パス、求解版SHA、入力SHA、
実行終了コード・確定時刻範囲・抽出法を併記して保持することを推奨する。
`complete=true` は実際の全計算完了時だけ指定し、時間切れやチェックポイント停止を変更しない。

|case|head_kind|coordinate|入力水頭|
|---|---|---|---|
|terzaghi|excess_pressure_head|down_from_drain|過剰水圧 / 9.81 kPa/m|
|moving_level|total_head|horizontal_from_left|全水頭|
|richards|total_head|up_from_water_table|H=ψ+z。ψだけを渡さない|

本体文書 `docs/RICHARDS.ja.md` の単位はs・m/s・m³/sである。
必要ならアダプター側で `time_h=time_s/3600`、`flux_m_h=flux_m_s*3600` と明示変換する。
積分流量Qを比流束qにするには実際の断面積・外向き符号を使う。
時間平均流束・工程lambda・セル平均水頭を瞬時の点値に偽装してはならない。
点値とセル平均は別比較であり、本契約は点値だけを対象とする。

### 質量収支

`balance` はメタデータJSONの任意配列。固定領域・源なし、断面積あたりの水量 m。

```json
"balance": [
  {"start": 0.02, "end": 0.2,
   "inflow_m": 0.0, "outflow_m": 0.003,
   "storage_start_m": 0.009, "storage_end_m": 0.006}
]
```

各量は本体の保存量と実際の受理刻みの流入出から抽出する。これは形式例。
流入・流出はともに非負の累積量、区間は連続、前区間末と次区間初めの貯留は
1e-12 m以内で連続とする（入力整合性検査、物理的受入許容値ではない）。
`defect_m=(S_end-S_start)-inflow+outflow` とその相対値を出す。
balance欠落は未評価。疎な瞬時流束CSVから水量を補間推定して閉じたとしない。
JOINTの領域/世代変更・撤去水量を含む契約は別途必要で、本版は対応しない。

### 水圧PDEと力学経路

任意 `mechanics` はTerzaghiの線形1Dだけ。各行は
`time, law, x, head, settlement_m`。`law` は
`linear_1d_mv_0.001_per_kPa_gamma_9.81`、`x` と `head` は全領域の節点列
（xの両端0/1、単調増加、3点以上）。沈下は圧縮正、m。
`--demo` がその例を生成する。
条件付き線形積分による経路誤差・水頭伝播誤差・標本積分誤差・合計を別出力する。
Richardsのこの欄への入力は拒否し、JOINTの独立3D力学が未評価であることを保つ。

本体CSV/JSONを抽出できた後の呼出し例（通常は約1秒〜数秒、点数依存）：

```powershell
& $python "$ref/compare.py" "$out/actual-head-flux.csv" `
  --metadata "$out/actual-metadata.json" --output "$out/actual-comparison.json"
if ($LASTEXITCODE -ne 2) { throw 'comparison did not produce expected undefined-acceptance diagnostic' }
```

合否基準は既存文書に定義された対象へだけ適用する。今回の過渡PDE誤差・質量収支・
力学経路の本体受入許容値は**未定義**なので、計測値が小さくても受入にしない。
A-04で必要な接合面水頭写像、世代/開閉、過渡水量、圧力荷重・仕事、結合時間細分は未達。
長時間3D・実CAD・一般Richards/JOINT求解は、本体入力と独立力学参照・受入基準が
整うまで実行手順も**未定義**。存在しない本体オプションや想定合格値を作らない。

## 既存Python全件と変更範囲の確認

既存テストとその閾値・分割器は変更していない。既存群には本体CLI・CAD依存があるので、
上記の独立最小環境とは別に、既存 `requirements-reference.lock.txt` とRust 1.99.1を用いる。
Linux実行結果の群別件数・ログは `results/validation.json` に記録する。
Windowsでは既存の次のコマンドを用いる（全件は多数の保存証拠再計測を含むため、
回線とSSD容量に依存して30〜120分以上の見込み、今回Windowsでは未測定）。
この長い全件確認と、数秒の新規独立参照を別々に計測する。

```powershell
py -3.13 -m venv work/py-water-regression-venv
if ($LASTEXITCODE -ne 0) { throw 'regression venv failed' }
$regPython = (Resolve-Path 'work/py-water-regression-venv/Scripts/python.exe').Path
& $regPython -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'regression dependency install failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'existing CLI build failed' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  & $regPython tools/run_ci_python_tests.py --group $group --output "$out/existing-$group"
  if ($LASTEXITCODE -ne 0) { throw "existing group failed: $group" }
}
git diff --name-status 2d2d54c48070581d088959b52984ceb348e0e85e HEAD
```

既存群は合計532件（各265/259/4/4）を収集する。
全群の `tests_run=assigned_total`、失敗・errors・skip・expected_failures・
unexpected_successesが全て0で初めて回帰成功。合計の欠落やskipで通さない。
Git差分はこの新規フォルダー内の `A` だけであることを確認する。
過去証拠の不足で失敗した場合は同じ起点の原ファイルを取得し、書換えず再実行する。
CIの時間切れ・キャンセルは成功に数えず、新先端の実際のrunと結論を別に確認する。

LinuxではPython 3.13.17と同じ既存依存lockで、次の4群を変更せず実行した。
実行環境で利用できるPythonを `python` とし、既存CLIを先にrelease buildする。

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
cargo build --locked --release
python tools/run_ci_python_tests.py --group python-core-0 --output work/python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/python-spatial-evidence
python interim/2026-11-19-py-transient-water-reference/checks.py -v
python interim/2026-11-19-py-transient-water-reference/diagnose.py --profile smoke --output work/smoke.json
```

各行の終了コードを別々に検査する。収集された532件と新規24件は別件数。
既存2証拠群は保存ファイルの再計測であり、新規の長時間NC/Richards/JOINT求解ではない。
