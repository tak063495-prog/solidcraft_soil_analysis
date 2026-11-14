# Windows PCでの再現と将来比較

作業はリポジトリーのルートからPowerShellで行う。出力は毎回新しいフォルダーへ保存する。参照生成・比較はRust/既存Pythonソルバを呼ばず、後半の「本体出力を取得する手順」だけが既存CLIを起動する。**P6-03/P6/NC-01は未達**。以下の所要時間はLinuxの小型参照約10秒以外は未実測の目安で、WindowsのCPU・メモリー・保存先に依存する。

## 依存と小型参照

既存の依存定義をそのまま使う。CIやlockを変更しない。Python 3.13を使用する。

```powershell
$Ref = 'interim/2026-11-21-py-mcc-independent-reference'
$Tag = Get-Date -Format 'yyyyMMdd-HHmmss'
$Run = "runs/py-mcc-independent-$Tag"
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
py -3.13 -m venv .venv-mcc
$Py = '.venv-mcc/Scripts/python.exe'
& $Py -m pip install --disable-pip-version-check -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'install failed' }
& $Py "$Ref/checks.py"
if ($LASTEXITCODE -ne 0) { throw 'independent software checks failed' }
& $Py "$Ref/run_reference.py" --output "$Run/reference"
if ($LASTEXITCODE -ne 0) { throw 'reference generation failed' }
& $Py "$Ref/assess_series.py" "$Run/reference/column-material-series-input.json" --output "$Run/material-diagnostic.json"
if ($LASTEXITCODE -ne 0) { throw 'material-series diagnosis rejected' }
& $Py "$Ref/assess_series.py" "$Run/reference/column-spatial-series-input.json" --output "$Run/spatial-diagnostic.json"
if ($LASTEXITCODE -ne 0) { throw 'spatial-series diagnosis rejected' }
```

目安：依存導入1〜10分、40新規試験1〜30秒、既定参照10秒〜2分、系列診断1秒程度。出力は全6材料点経路の128/256/512生系列・閉形式等方圧密・元MC柱・特定したMCC柱・交換JSON/CSV・出典hash・`summary.json`。正常終了0は生成成功であり数値受入ではない。比較許容値は未定義、`accepted=null`。

MCC柱の材料系列は256/512/1024圧力増分、GL16固定。深さ系列は8/16/32セル、材料2048固定。GL16→32は同じ材料2048で測定する。全体施工刻み・材料内分割・空間細分を同時に変えた系列から、空間次数を認定しない。

さらに細かい**独立1Dだけ**の確認は別出力で実行できる（目安1〜20分、未実測）。

```powershell
& $Py "$Ref/run_reference.py" --output "$Run/finer-reference" `
  --levels 512 1024 2048 `
  --column-material-levels 1024 2048 4096 `
  --column-spatial-cells 32 64 128 --column-spatial-material-n 8192
if ($LASTEXITCODE -ne 0) { throw 'finer reference failed' }
```

判定：生の差、観測次数、外挿の条件、降伏/NCL/体積整合ドリフトを確認する。次数がない定数・反転・丸め域は「未定義/未確認」。材料点の終点をS2の漸近限界値と直接比較して合格にしない。GLやRichardsonを厳密解と呼ばない。新規の数値許容値は作らず、未定義ならnullを保持する。

## 将来の本体CSV/JSONとの比較

本体側で必要な全履歴・初期状態・単位・測度が出力できてから実行する。参照JSONの`metadata`と同じ物理入力・点/節点・負荷座標を明示する。`source_sha256`は本体実装/実行版の実SHAに置き換え、`physics_sha256`は同じ入力定義を検証したうえで一致させる。値を合わせるために参照メタデータを無検証でコピーしない。

CSV数値は17桁の有効数字を推奨する。表示用に丸めた値を完全精度の履歴としない。

最小CSV列：

```
GP:   step,id,load_coordinate,p,q,pc,e,eps_v
Node: step,id,load_coordinate,ux,uy,uz
```

JSONは`{"metadata": {...}, "rows": [...]}`。CSVは同じメタデータだけを別JSONに保存する。`schema=independent-mcc-history-v1`、kPa/m/無次元、有効応力圧縮正、εvの圧縮正と`strain_measure=log_volume`、変位はglobal_xyzを明示する。生成例は`k0-ocr1-gp-reference.*`、下層だけの`documented-mcc-column-gp-reference.json`と`documented-mcc-column-node-reference.json`。

```powershell
# 以下の本体ファイルは将来の明示出力。未出力の列を推測して作らない。
$BodyGp = 'runs/my-body-export/gp.csv'
$BodyGpMeta = 'runs/my-body-export/gp.metadata.json'
& $Py "$Ref/compare_history.py" `
  --reference "$Run/reference/documented-mcc-column-gp-reference.json" `
  --candidate $BodyGp --candidate-metadata $BodyGpMeta `
  --output "$Run/compare-column-gp.json"
if ($LASTEXITCODE -ne 0) { throw 'GP comparison rejected' }
& $Py "$Ref/compare_history.py" `
  --reference "$Run/reference/documented-mcc-column-node-reference.json" `
  --candidate 'runs/my-body-export/nodes.json' `
  --output "$Run/compare-column-nodes.json"
if ($LASTEXITCODE -ne 0) { throw 'node comparison rejected' }
```

目安：小型出力1秒〜1分、100万行はメモリーやJSONサイズに依存し数分以上。現比較器は全データをメモリーへ読むため、大容量はPCで範囲を明示した完全なID集合ごとのファイルに分ける。途中の履歴行を間引いたものを「全履歴」としない。

出力は各量の最大絶対差・相対L2・元ノルム・入力hash。ゼロ参照ノルムの相対差はnull。数値許容値未定義、全てaccepted=null。拒否は終了2、出力を生成しない。既存ファイルへ上書きしない。

pcが減少するOCR4膨張枝は、モデル結果を保持したまま今回の単調pc交換契約の対応外として拒否する。これはその材料経路を物理的に不正と判定する意味ではない。p'<=0、単位/符号/有効・全応力の不一致、欠測・未定義履歴の補完、NaN、pc逆行、e/εv不整合、不完了runは拒否。工程終点しか出せない本体出力を内部全履歴へ補完しない。真のεvが出ていなければ、材料塑性ひずみや変位から勝手に置き換えない。

## 3D本体出力を取得する手順（PCのみ）

以下はS8/S9の**乾燥MCC柱**だけで、C-03や湿潤JOINT等の並行作業は起動しない。生成される本体結果の最大32層は192 Tet10だが、3段階×複数積分制御の独立した空間系列であり、クラウドでは再求解しない。既存の大規模CAD壁・10万要素系列も今回のクラウド実行対象にしない。

まず既存の入力準備器だけを呼び、同じ物理条件の2/4層入力を作る。これは本体比較の入力準備であり、新規独立参照の生成経路ではない。

```powershell
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'body CLI build failed' }
$Solver = (Resolve-Path 'target/release/soil-fem-3d.exe').Path
$Inputs = "$Run/body-inputs"
& $Py tools/prepare_mixed_mcc_reference.py $Inputs
if ($LASTEXITCODE -ne 0) { throw 'input preparation failed' }
# 既存の幾何作成関数だけを使い、S9と同じ物理を8/16/32層へ展開する。
@'
import json,sys
from pathlib import Path
sys.path.insert(0,'tools')
from prepare_mixed_mcc_reference import column
root=Path(sys.argv[1])
for n in (8,16,32):
    d=root/f'n{n}';d.mkdir(exist_ok=False)
    (d/'mesh.json').write_text(json.dumps(column(n)),encoding='utf-8')
    (d/'stages.json').write_bytes((root/'n4/stages.json').read_bytes())
# 空間軸を固定刻みで3水準。別に、16層固定の全体/材料刻み対照を3水準ずつ。
for label,n,load,material in (
    ('space8',8,1/32,1e-5),('space16',16,1/32,1e-5),('space32',32,1/32,1e-5),
    ('load16',16,1/16,1e-5),('load32',16,1/32,1e-5),('load64',16,1/64,1e-5),
    ('mat1e4',16,1/32,1e-4),('mat1e5',16,1/32,1e-5),('mat1e6',16,1/32,1e-6)):
    d=root/label;d.mkdir(exist_ok=False)
    (d/'mesh.json').write_bytes((root/f'n{n}/mesh.json').read_bytes())
    case=json.loads((root/f'n{n}/stages.json').read_text())
    case['settings'].update(initial_increment=load,material_maximum_substep=material)
    (d/'stages.json').write_text(json.dumps(case),encoding='utf-8')
'@ | & $Py - $Inputs
if ($LASTEXITCODE -ne 0) { throw 'series preparation failed' }
foreach ($Case in @('space8','space16','space32','load16','load32','load64','mat1e4','mat1e5','mat1e6')) {
  & $Solver solve-stages "$Inputs/$Case/mesh.json" "$Inputs/ground.json" `
    "$Inputs/$Case/stages.json" "$Inputs/$Case/result.json" 2>&1 | Tee-Object "$Inputs/$Case/solve.log"
  if ($LASTEXITCODE -ne 0) { throw "body solve failed: $Case (not collapse/pass)" }
}
```

目安：release build2〜15分、9経路は30分〜数時間、未実測。疎反復上限・残差未達・切戻し枯渇はそのまま停止してログを保存する。最大反復数、残差基準、pc/eの初期化、形状、K0、材料、載荷・撤去は調整して通さない。caseの`initial_increment`は実際の受理増分を保証しないので、必ず受理列と材料内分割を結果から確認する。自動cutbackで意図した系列にならなければ「未確認」。

この手順の現本体JSONに将来比較に必要な全GPのp/q/pc/e/εvと共通物理点がすべて揃うとは限らない。特に非一致メッシュの共通物理点写像、全体積ひずみ測度、応力6成分の奥行き比較の本体exportは**未確認/別途必要**。既存Pythonソルバで不足履歴を作って本体出力にしない。元S8の当時のpc不合格を今回の別版・別出力で上書きしない。

押出しMC柱の3D h/奥行き対照も必要ならPCで既存入口を使う。これはMCCと別物理で、元5%/10%/1%判定を保持する。

```powershell
& $Py tools/plastic_spatial_convergence.py "$Run/mc-inplane" --binary $Solver `
  --divisions 4 8 16 --depth-divisions 4 --increment 0.03125 --material-substep 0.00001 --timeout 7200
if ($LASTEXITCODE -ne 0) { throw 'MC h series incomplete' }
& $Py tools/plastic_spatial_convergence.py "$Run/mc-depth" --binary $Solver `
  --axis depth --divisions 4 8 16 --in-plane-divisions 8 `
  --finite-width --increment 0.03125 --material-substep 0.00001 --timeout 7200
if ($LASTEXITCODE -ne 0) { throw 'MC depth series incomplete' }
```

目安は各数十分〜数時間、未実測。ここでの既存MC掘削形状は元.2 m一様柱の再実行ではないため、6.828%/1.034%の書換えやMCC空間次数には使わない。Griffithsや湿潤JOINTの系列へ流用しない。有限幅結果・深さ系列が面内h系列の代用にはならない。

## 将来の系列診断と合否

`column-*-series-input.json`が入力例。ルートに`axis`（spatial/load_increment/material_increment）、`parameter_unit`（spatialはm、増分は明示した1又はkPa）、`run_status=completed`。3つの`levels`には物理入力hash、量/単位、共通sample_ids、値、深さ、拘束、負荷経路、全体荷重/材料刻み、mesh_sha256を保存する。空間軸は`h_m`と`nested_h_refinement=true`も必須。

本体の独立系列を同じ契約へ**補間なしで**明示出力できてから、次を実行する。

```powershell
& $Py "$Ref/assess_series.py" 'runs/my-body-export/spatial-series.json' --output "$Run/body-space-diagnostic.json"
& $Py "$Ref/assess_series.py" 'runs/my-body-export/load-series.json' --output "$Run/body-load-diagnostic.json"
& $Py "$Ref/assess_series.py" 'runs/my-body-export/material-series.json' --output "$Run/body-material-diagnostic.json"
```

3水準・同じ物理点でないデータ、深さ/拘束/荷重/材料を同時変更した空間系列は拒否する。原5%等が適用される限定試験だけは原定義で別判定する。新しいMCC材料点/1D参照や押出し・有限幅応力への数値許容値は既存文書で未定義なのでaccepted=null。失敗・判定不能・未達・キャンセル・時間切れを合格にしない。

## Linux/Windowsの既存Python全件確認

既存CIと同じ4群を使う。NC証拠各4試験は保存済み大きなarchiveを読むので、該当partを全て取得しておく。旧証拠の再監査であり大規模MCCの新規求解ではない。`interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part*`だけは現CIと同じ除外が可能だが、mean/spatialのpartを除外しない。

```powershell
& $Py "$Ref/validate_linux.py" --output "$Run/validation"
if ($LASTEXITCODE -ne 0) { throw 'Python validation failed or incomplete' }
```

この入口はWindowsでも使える。名前のlinuxは保存した最終実測環境を示す。目安10〜60分以上、証拠の復元を含む読込で一時空間10〜30 GBを見込む（未実測目安）。`validation.json`の各群tests_run、fail/error/skip/expected-failure/unexpected-success=0、全発見IDの一回実行、前後ソースhash一致、new_checksの40件を確認する。途中停止ならstatusがrunning/failedのままで、成功とはしない。CI設定・発見ID・件数条件を変更しない。
