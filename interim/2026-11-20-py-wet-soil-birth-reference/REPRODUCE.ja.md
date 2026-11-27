# Windows PCでの再現と将来の本体照合

PowerShell 7、Python 3.13、Git、既存lockのNumPy/SciPy環境を使う。
以下はリポジトリルートから実行する。作業前に未保存の変更を保存する。
新ブランチは指定SHAの直接の子で、PR #13の後続更新や並行draft PRを混ぜない。

```powershell
git fetch origin codex/py-wet-soil-birth-reference
git switch codex/py-wet-soil-birth-reference
git rev-parse HEAD
git log -1 --format=%P
# 親に 2d2d54c48070581d088959b52984ceb348e0e85e を確認

if (!(Test-Path .venv\Scripts\python.exe)) { py -3.13 -m venv .venv }
$py = (Resolve-Path .venv\Scripts\python.exe).Path
& $py -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw '依存導入失敗' }
$env:OPENBLAS_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:PYTHONUTF8='1'
$env:PYTHONIOENCODING='utf-8'
$ref = 'interim/2026-11-20-py-wet-soil-birth-reference'
$run = 'runs/NEW-wet-soil-birth-reference-' + (Get-Date -Format yyyyMMdd-HHmmss)
New-Item -ItemType Directory -Force $run | Out-Null
& $py --version
& $py -m pip freeze | Out-File "$run/packages.txt"
```

## 小型参照と明示試験

```powershell
& $py "$ref/checks.py" 2>&1 | Tee-Object "$run/new-checks.log"
if ($LASTEXITCODE -ne 0) { throw '独立試験失敗' }
& $py "$ref/run_reference.py" --output "$run/reference"
if ($LASTEXITCODE -ne 0) { throw '参照生成失敗' }
```

目安：試験1〜10秒、参照10〜60秒（PC性能依存、実測を保存）。Linuxの小型参照は約8.5秒。
出力は7ケースのprofiles.json、全節点力・面圧・水頭・流束JSON、measurements.json、synthetic-onlyのJSON/CSV束。
検査条件は`checks.py`に固定。失敗・エラー・skip・xfail・unexpected successは0であること。
追加44件を全件実行する（既存532件と合わせて576件）。
参照のexit 0は生成完了だけを意味し、A-03の合格ではない。
丸め誤差域の次数はnull、非適合水頭系列は約1次、生の切断要素流束の最大差は収束未確認として残す。
生成JSONの`accepted=null`、`a03_complete=false`を維持する。

## 既存Python全件（必要証拠の完全取得後）

既存CIはpytestではなくunittest。新試験は既存収集に入らない`checks.py`へ置いている。
通常のcloneでは既存のZIP/分割アーカイブも必要。sparse checkoutを使う場合は、少なくとも既存試験の参照先を含める。
封印NCのmean-stress約5.49GB・spatial約0.92GBと、湿潤JOINT等の証拠を使う。
不足ならGitから**同じ起点のblob**を取得する。参照値・SHA・失敗記録を作り直して埋めない。
空き容量は既存履歴の展開を含め30GB以上を目安にする。
`interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*`は既存CI同様、今回4群に不要。

```powershell
# 既存Python試験が使うCLI。独立参照器自身はこのexeを呼ばない。
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw '既存CLIのbuild失敗' }
& $py "$ref/validate.py" --output "$run/validation" 2>&1 | Tee-Object "$run/validate.log"
if ($LASTEXITCODE -ne 0) { throw '全件検証失敗。ログを保持する' }
Get-Content "$run/validation/validation.json"
```

目安：build数分、全件10〜30分以上。封印証拠の取得時間は回線依存。
大規模な新規求解ではなく、既存全件試験と既存封印データの再計測である。
各群のログ/完全ID、追加試験ID、前後のソースSHAを保存する。
既存532件（core-0:265、core-1:259、mean:4、spatial:4）と追加試験がすべて実行され、
fail/error/skip/expected failure/unexpected success=0、発見IDの欠落・重複0、source_hashes_unchanged=trueが条件。
NCの封印証拠試験成功をNC-01の数値受入へ転用しない。
時間切れ・キャンセル・途中ログは未完了であり成功件数へ入れない。

## PC用の長めの独立系列（クラウドでは未実行）

```powershell
& $py "$ref/run_reference.py" --pc-series --output "$run/reference-long"
if ($LASTEXITCODE -ne 0) { throw '長系列未完了' }
```

16/32/64/128/256/512を、適合・非適合、上向き・下向きで別々に実行。
目安1〜10分、出力measurements.jsonと同じケース別JSON。線形sparse solve自体は小さいが、全節点力の分割積分も各系列で行う。
空間と積分次数を混ぜない。各誤差・観測次数、水頭境界流入出、全節点力残差、合力・モーメントを比較する。
過去の小型参照値は上書きせず、新しい`runs/NEW-*`へ保存する。
本体用の新しい許容値や合格数値を系列から逆算しない。全体acceptedはnullのまま。

## JSON / CSV比較器の動作確認

まず付属の**人工例**を往復する。これは本体比較ではない。

```powershell
& $py "$ref/compare.py" --case partial_new_layer `
  --input "$run/reference/synthetic-only/export.json" `
  --history-before "$run/reference/synthetic-only/history-before.json" `
  --output "$run/synthetic-json-comparison.json"
if ($LASTEXITCODE -ne 2) { throw '未定義許容値の正常診断はexit 2' }

& $py "$ref/compare.py" --case partial_new_layer `
  --input "$run/reference/synthetic-only/csv" `
  --history-before "$run/reference/synthetic-only/history-before.json" `
  --output "$run/synthetic-csv-comparison.json"
if ($LASTEXITCODE -ne 2) { throw 'CSV診断失敗' }
```

出力の各誤差は0、ただしvalid_export=trueは形式・出典検査であり、accepted=null。
数値差が大きくても許容値未定義のためexit 2の診断を出し、合格にはしない。
形式・由来の拒否はexit 1、accepted=false。

## 将来の本体出力アダプター

現時点で本体A-03拡張は範囲外。成功する未定義のBIRTH CLIを創作しない。
既存の静水圧MC受入入力をDarcy/MCC/部分水没の受入に流用しない。
本体に明示仕様・対応出力が追加された後、ユーザーPCで各小型ケースを実行し、以下の契約へ変換する。

JSONはsynthetic-only/export.jsonを**構造の例**に使う。
CSV束はmetadata.json、materials.json、nodes.csv、forces.csv、fluxes.csvの5ファイル。

|項目|必須契約|
|---|---|
|metadata|schema、case、units、head_kind=total_z_up、stress_sign=compression_positive、complete=true、status=completed、generation=1、stage=birth_initial_state、4種の補完/転用フラグ=false|
|nodes|id,x,z,body,head_applicable,head_m,pressure_kpa,effective_vertical_kpa,effective_horizontal_kpa,pressure_origin。ID・座標は付属Q4参照と一致。非適用水頭はnull/空欄|
|forces|各節点のpore_internal/effective_internal/body/pressure_exterior/pressure_interface/effective_boundary/gradient_volume/seepageの各_x,_z。kN、単位奥行|
|fluxes|boundary=bottom/top、外向き正outward_m3_s。節点反力を面境界別に合計し、元の節点流量も別保存|
|materials|oldの全状態、新設全点のid,model,origin=explicit_fresh_input,source_id=null,plastic_strain[6],equivalent_plastic_strain,plastic_work_kpa,last_return。MCCはeffective_stress_kpa[6],pc_kpa,void_ratio,Mも明示|
|history-before|出力生成前に独立保存したold全状態とnew_ids全ID。出力から逆生成しない|

参照の節点IDは旧土0〜9、新土10〜19（各層4 Q4、x=0/1）。両側の境界節点は別ID。
Tet10/3Dなど異なるメッシュを無断補間して同じ参照へ比較しない。
アダプターはこの矩形の整合節点力を別途明示的に投影し、その仮想仕事と投影誤差を検証する必要がある。
現比較器は同一参照メッシュの入力だけに対応し、任意メッシュの受入を主張しない。
水圧を節点へ無断外挿・補完して契約を満たしてはならない。GP水圧しかない場合は明示出力が揃うまで拒否する。
material snapshotの点集合は入力manifestで全点列挙し、人工例の1点を実要素の全GPと見なさない。
PRESCRIBED内部貯水槽例の比較は現CSV契約の上下端2境界だけでは足りないため、本体比較未対応。
その場合はreservoir_outを含む拡張契約を先に仕様化する。値を上下端へ紛れ込ませない。

```powershell
# 下記入力は将来の明示アダプター出力。今は未生成。
& $py "$ref/compare.py" --case darcy_upward `
  --input runs/NEW-production-export/export.json `
  --history-before runs/NEW-production-input/history-before.json `
  --output "$run/production-comparison.json"
# exit 1=拒否、exit 2=診断完了・accepted=null。exit 0を成功条件にしない。
```

本体の実計算・適合メッシュ・全材料点出力・保持力/水頭世代・再開の仕様が揃うまでA-03未達。
新設MCC・古い履歴進展の物理仕様と適用許容値は別途定義・レビューが必要。
Windowsの新先端CI成功はここでは未確認。Linux成功・旧CI・時間切れ/キャンセルを転用しない。
