# Windows PCでの再現

PowerShell 7、Git、Python 3.13。既存作業ツリーを切り替えず、新しいcloneを使う。
時間はPC向けの未実測の計画値。出力先は毎回新規にし、途中・失敗を上書きしない。

```powershell
git clone --single-branch --branch codex/py-wet-soil-birth-darcy-reference https://github.com/tak063495-prog/solidcraft_soil_analysis.git C:\soil-darcy-ref
if ($LASTEXITCODE -ne 0) { throw 'clone failed' }
Set-Location C:\soil-darcy-ref
git merge-base --is-ancestor 702a19a7ce189f1173caed727ee2a0059298a607 HEAD
if ($LASTEXITCODE -ne 0) { throw 'wrong base' }
git rev-parse HEAD
git log -1 --format=%P
# このPRの親が702a19a7ce189f1173caed727ee2a0059298a607であることを確認
py -3.13 -m venv .venv-ref
if ($LASTEXITCODE -ne 0) { throw 'venv failed' }
$py = (Resolve-Path .venv-ref\Scripts\python.exe).Path
& $py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'install failed' }
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$ref = 'interim/2026-11-26-py-wet-soil-birth-darcy-reference'
$run = 'runs/darcy-birth-' + (Get-Date -Format yyyyMMdd-HHmmss)
if (Test-Path $run) { throw 'output must be new' }
New-Item -ItemType Directory $run | Out-Null
git rev-parse HEAD | Out-File "$run/head.txt"
& $py -m pip freeze | Out-File "$run/packages.txt"
& $py "$ref/checks.py" 2>&1 | Tee-Object "$run/new-tests.log"
if ($LASTEXITCODE -ne 0) { throw 'new tests failed' }
& $py "$ref/reference.py" --output "$run/reference"
if ($LASTEXITCODE -ne 0) { throw 'reference incomplete' }
```

小型参照・新規42試験は数秒〜1分を計画。出力はα2種×n=8/16/32/64の8 JSONとconvergence.json、policy-audit.json。
各JSONは全水理点、全機械点、全節点の力成分、境界流量、誤差・残差・位相監査を含む。
生成・試験成功とA-03受入を区別する。計算の健全性試験はfail/error/skip/xfail=0。
水位/水頭の約2次、収支、残差の全生値を確認するが、本体受入許容値は未定義、accepted=nullを維持する。
根求解丸め差の水面負圧を0へ書換えない。低次数・非単調・未完了はそのまま保存する。

## JSON/CSV契約の人工往復

```powershell
$example = "$run/reference/alpha0.5-n8.json"
& $py "$ref/compare.py" --input $example --reference $example --output "$run/synthetic-json-diagnostic.json"
if ($LASTEXITCODE -ne 2) { throw 'diagnostic must finish with undefined acceptance (2)' }
& $py -c "import sys;sys.path.insert(0,sys.argv[1]);from compare import read_json,write_csv;write_csv(read_json(sys.argv[2]),sys.argv[3])" $ref $example "$run/synthetic-csv"
if ($LASTEXITCODE -ne 0) { throw 'CSV generation failed' }
& $py "$ref/compare.py" --input "$run/synthetic-csv" --reference $example --output "$run/synthetic-csv-diagnostic.json"
if ($LASTEXITCODE -ne 2) { throw 'CSV diagnostic failed' }
```

所要1分以内、JSON/CSVの全差0だがaccepted=null。本体比較ではない。
CSV束はmetadata.json＋hydraulic.csv/samples.csv/forces.csv/budget.csv。
単位・符号・出典・状態・ID・座標はREADMEの契約に従う。
明示的な本体アダプターが未実装であるため、本体の未対応CLI・入力schemaを創作しない。

## 既存Python全件

独立参照にRustは不要だが、既存試験はRust CLI/examplesを使用する。
通常の完全cloneに既存ZIPとNC分割証拠が必要。欠けていれば同じ起点のGit blobを取得する。
保存証拠の再計測であり、新しいNC大型求解ではない。空き容量30 GB以上、通信量は数GB以上を見込む。

```powershell
$env:CARGO_BUILD_JOBS='2'
$env:RUST_TEST_THREADS='2'
cargo build --locked --release --bins --examples
if ($LASTEXITCODE -ne 0) { throw 'existing CLI build failed' }
$groups = 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'
$codes = @{}
foreach ($g in $groups) {
  & $py tools/run_ci_python_tests.py --group $g --output "$run/$g" 2>&1 | Tee-Object "$run/$g-process.log"
  $codes[$g]=$LASTEXITCODE
}
$codes | ConvertTo-Json | Out-File "$run/existing-exit-codes.json"
# 全4群の実行後に判定する。途中の失敗で残り群を省略しない。
if (@($codes.Values | Where-Object { $_ -ne 0 }).Count -gt 0) { throw 'existing suite has errors; preserve all logs' }
```

build5〜20分、全件10〜60分＋証拠取得時間を計画。既存568 IDを4群276/284/4/4へ一意割当。
件数が変わればID一覧を比較し、収集条件を書き換えて合わせない。
旧NCのLinux setUpClass準備失敗（meanのsaved mean-stress/complement partition差、spatialのremeasured JSON差）は
各4メソッド未実行、成功にもskipにも数えない。Windowsでの成功は新しい実行でのみ確認する。
Linuxだけplatform許容を追加したり旧値・厳密比較を変更したりしない。
全件成功の条件は全568本体実行・fail/error/skip/expected-failure/unexpected-success=0・全群exit0。
一部成功や時間切れ・キャンセルは全件成功ではない。

## PC専用の長い1D系列

```powershell
& $py "$ref/reference.py" --pc-series --output "$run/reference-long"
if ($LASTEXITCODE -ne 0) { throw 'long series incomplete' }
```

α2種×n=32/64/128/256/512/1024。計画1〜10分、同形式の12 JSONとconvergence.json。
クラウドではこの系列を実行しない。小型値を上書きせず、誤差の丸め飽和・次数劣化も保持する。
本体受入値を最細結果や外挿から逆算しない。accepted=nullのまま。

## 将来の本体比較と大きな3D系列

本体側がA-03用初期水理・新設MCC/旧履歴・全出力を明示仕様化した後にだけ、
元入力・元出力・実行版SHAを固定した別runで実行し、明示アダプターから同一IDのJSON/CSVを作る。
生成参照を編集したり、水圧/節点力を黙示補間したりしない。

```powershell
# 下のinputは将来の明示アダプター出力。現時点では未生成。
& $py "$ref/compare.py" --input runs/NEW-production-export --reference "$run/reference/alpha0.5-n32.json" --output "$run/production-diagnostic.json"
if ($LASTEXITCODE -eq 1) { throw 'incompatible/missing export; do not interpolate' }
if ($LASTEXITCODE -ne 2) { throw 'unexpected diagnostic exit' }
```

想定1分以内。差・節点残差・水収支を出すが適用許容値未定義のためaccepted=null。
本体の水圧や応力は全GP出力と本参照節点が異なる場合、このコマンドでは比較できない。

大きな3Dは**この参照器では未実装**。未存在の3Dコマンドを成功手順として載せない。
将来、明示仕様と独立3D参照・同一ID出力がそろった時点でPC上の3水準を固定し、
少なくとも水頭/自由面/応力/全節点力/各境界流量/全材料履歴、平衡・質量・再開を保存する。
計画30分〜数時間/水準、メモリー・実測時間は未確認。
3D実行ファイル・固定メッシュ・初期履歴・比較許容値が今は無いため、実行可能な3D受入一式は未達。
NC-01/P6/A-03の完了、既存限定受入の拡張、新先端Windows CI成功として扱わない。
