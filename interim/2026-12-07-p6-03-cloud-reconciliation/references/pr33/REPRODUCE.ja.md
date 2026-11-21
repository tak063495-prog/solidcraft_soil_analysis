# Windows PC での再現

指定起点 `3b1133d3be9905dc24f62c27b5b1eb456a83f2e2` から分岐したこのブランチを使う。保存先の日付は依頼ラベルで、実施日は2026-11-12（日本時間）。マージ不要。以下はリポジトリのルートで PowerShell を実行する。既存の出力先は再利用しない。

```powershell
git switch codex/py-p6-03-extrusion-vibration-depth-reference
$ref = 'interim/2026-11-28-py-p6-03-extrusion-vibration-depth-reference'
py -3.13 -m venv .venv-p6-reference
$py = '.venv-p6-reference/Scripts/python.exe'
& $py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency installation failed' }
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYTHONUTF8 = '1'
New-Item -ItemType Directory -Force runs | Out-Null
& $py "$ref/collect_legacy.py" runs/p6-03-legacy-pc-v1.json
if ($LASTEXITCODE -ne 0) { throw 'legacy evidence differs; do not rewrite it' }
& $py -m unittest discover -s $ref -p test_reference.py -v
if ($LASTEXITCODE -ne 0) { throw 'new reference tests failed' }
& $py "$ref/run_reference.py" runs/p6-03-small-pc-v1
if ($LASTEXITCODE -ne 0) { throw 'reference series incomplete' }
```

`runs/` が無い新規配置では先に `New-Item -ItemType Directory runs`。新規参照には Rust 実行版も既存 Python ソルバも不要。固定依存の NumPy/SciPy と標準ライブラリだけを使う。新しい試験はこのフォルダー内に置き、既存 `tests/`・`tools/` の収集対象・workflow・件数条件を変更しない。

小型系列は一般的な近年のPCで約1〜5分、メモリー目安1 GB未満（保証値ではない）。Linux実測は `results/summary.json`。Lamé は各2D/3D・平面ひずみ/GPSの3水準、振動は5周期・20/40/80/160刻み、有限箱は幅1/2 m×奥行き2/4/8 m×メッシュ2水準。出力は `summary.json` とケース別 JSON/CSV/メタデータ JSON。物理点・単位・成分を各レコードに保持する。

`summary.json` の `accepted=null`・`p6_03_complete=false` が正しい。元の比0.7・次数範囲・符号条件は旧MC系列の参照なし外挿の必要条件であり、新しい弾性・動解析の受入許容値ではない。新規ベンチマークの本体比較許容値は **未定義**。unit test の丸め誤差や既知次数の検査をP6の合否基準へ使わない。旧3項目は元の判定不能のまま。

## 大きな系列と長い応答（ユーザーPC専用、クラウド未実行）

小型系列の後、空きメモリーを確認して一つずつ実行する。下記3Dは直接疎行列分解なので、メッシュ拡大に対してメモリーは要素数比例では増えない。

```powershell
# まず単独ケースで費用を確認：約数分〜30分、数GBの可能性
& $py "$ref/run_reference.py" runs/p6-03-pc-mesh-v1 --widths 1 --lengths 4 --levels 1 2 4 --ring-levels 4 8 16
if ($LASTEXITCODE -ne 0) { throw 'incomplete; no acceptance or extrapolation' }
# 奥行き増加を固定メッシュごとに観測：約5〜60分、数GB〜16GB以上の可能性
& $py "$ref/run_reference.py" runs/p6-03-pc-length-v1 --widths 1 --lengths 4 8 16 32 --levels 1 2 --ring-levels 4 8 16
if ($LASTEXITCODE -ne 0) { throw 'incomplete; no plane-strain-limit acceptance' }
# 3D費用を増やさず100周期の位相誤差・散逸を検査：約1〜10分
& $py "$ref/run_reference.py" runs/p6-03-pc-long-v1 --widths 1 --lengths 2 --levels 1 --ring-levels 4 8 16 --periods 100
if ($LASTEXITCODE -ne 0) { throw 'long response incomplete' }
```

上の所要時間は未実測の概算で、上限や性能保証ではない。メモリー不足・中断時は不完了として出力先を保持し、次は新しい出力先を使う。最後まで到達した `summary.json` だけを完了系列とする。部分CSV・時間切れ・キャンセルを成功として集計しない。

メッシュ3水準以上でも有限箱の荷重端の特異性、固定底面と側面の影響が残る。中央の同一物理点の差と、同じメッシュの2D平面ひずみとの差、自由端と中央の差を分ける。`mesh_differences` は差のみで外挿しない。追加水準の次数を調べる場合も同じB/L・同じ物理点・同じ片側トレースを使い、長さ変更とh変更を混ぜない。平面ひずみ極限の一般精度ゲートは未定義。

## 将来の本体 CSV/JSON との比較

本体の exporter/schema は追加していない。将来の出力から、単位・軸・全/有効応力・符号・成分・物理条件・点ID・位置・物理時刻を出力仕様に照らして明示的に対応付けた交換ファイルを**別途**用意する。本体計算もその対応付けも今回未実行。旧MCの y=奥行き・z=鉛直・圧縮正kPaを、新参照の z=奥行き・引張正Paへ黙って読み替えない。時刻0の静的状態と物理0秒も区別する。

JSONの例は `lame-plane-strain.json`、時刻歴は `sdof-history.json` / `rod-history.json`。CSVは同名の `.meta.json` を明示指定する。正しい単位は応力 Pa、変位 m、速度 m/s、加速度 m/s^2、周期 s、振動数 Hz、エネルギー J。JSON/CSVの完全な列と条件は `compare_outputs.py` の `FIELDS` / `validate` および実際の出力で確認できる。CSVはその列順、UTF-8、各列必須。

```powershell
# 交換形式の自己往復の例。本体比較の証拠ではない。
& $py "$ref/compare_outputs.py" runs/p6-03-small-pc-v1/sdof-history.json runs/p6-03-small-pc-v1/sdof-history.csv runs/p6-03-roundtrip-v1.json --candidate-manifest runs/p6-03-small-pc-v1/sdof-history.meta.json
if ($LASTEXITCODE -ne 2) { throw 'expected diagnostic-only exit 2' }
# 将来用。candidate.json は本体の根拠付き対応付けを完了してから用意する。
& $py "$ref/compare_outputs.py" runs/p6-03-small-pc-v1/lame-plane-strain.json runs/body-lame-candidate.json runs/p6-03-body-comparison-v1.json
if ($LASTEXITCODE -eq 1) { throw 'exchange invalid; no comparison verdict' }
if ($LASTEXITCODE -eq 2) { Write-Host '差だけを保存。許容値未定義、accepted=null' }
```

比較器の終了コード1は入力拒否、2は差の診断のみ。0を合格として返す経路はない。欠落/重複ID・時刻・成分、余分な点、座標差、単位差、符号/軸/応力基準の宣言差、物理条件差、非有限値、不完了を拒否する。補間・単位変換・符号反転・不足値補完は行わない。同じ宣言のまま値だけが逆符号/入替された場合は差を保存する。入力の宣言が真実かを一般的に証明する機能はなく、元の出力仕様の確認が必要。

## 既存Python全件を別途実行

既存試験は本体を使う。独立参照の生成経路とは分け、元の収集器を使って4群をすべて実行する。事前に参照の結果を上書きしない。

```powershell
cargo build --locked --release --bins --examples
if ($LASTEXITCODE -ne 0) { throw 'existing executable build failed' }
foreach ($group in @('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')) {
  & $py tools/run_ci_python_tests.py --group $group --output "runs/p6-03-existing-pc-v1/$group"
  Write-Host "group=$group exit=$LASTEXITCODE"
}
```

全群を終了コードとともに記録する。Linuxで既知の旧NC2群には厳密比較の不一致があるが、**今回の実測は分割証拠未配置によるsetUpClassのFileNotFoundError**であり、厳密数値比較を再実行した証拠ではない。各4メソッド未実行を成功・skipへ変更しない。完全証拠（mean約5.491 GB、spatial約.920 GB）を元コミットのblob/SHAと照合して配置してから、その厳密比較を再評価する必要がある。旧値・SHA・ソース・platform規則を直さない。Windowsの実結果はWindowsで確認し、Linuxや過去CIの結果を転用しない。今回の新先端CIは `results/validation.json` の記載以上に成功扱いしない。
