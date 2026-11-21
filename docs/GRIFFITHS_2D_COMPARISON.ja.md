# Griffiths斜面の661節点実データと5％比較

2026-11-08。起点はPR #13の `9187211d4fa224d7eaf936bf73818b06fb9fe190`。生産CLIの数値経路、元の閾値と過去の不合格記録は変更していない。クラウドでは1200／4800／19200 Tet10の求解を実行していない。新しい面内h系列の求解と5％達成は未完了であり、一般3D空間収束・P6全体も未完了である。

## 実データの所在と出典

リポジトリ全体の検索で、元の変位は **存在する** と分かった。`docs/inconsistent-evidence24` 自体はこのリポジトリに無いが、次の既存JSONに `source2d_u_xz_m` の661点×2成分が残っている。

- `interim/2026-10-30-Root635-full-audits-slope-reference/slope/coarse-field-samples.json`
- `interim/2026-10-30-Root635-CI-and-refinement-controls/slope/{xy1-ny1,xy1-ny2,xy1-ny1-global64}/all661-nodes-and200-prisms.json`
- `interim/2026-10-30-async-UI-and-independent-SRM-bound/slope/material-half1200/all661-nodes-and200-prisms.json`

起点の追跡2639ファイル、`docs`、`interim`、`tests/data`を検索し、101 ZIPの30058メンバーも名前とJSON/Python/Markdown/テキスト/CSV内容を検索した。元の生結果ファイルは直接同梱されていなかった。上の5コピーの変位は全て一致し、さらに固定版 `tak063495-prog/excel-vba-soil-fem@440d246d9e32f5892c79539d4a7b797ad29442ee` の実アーカイブを取得して全点の完全一致を確認した。Excelは再実行していない。

|出典|SHA-256（元ファイルのバイト列）|
|---|---|
|`docs/inconsistent-evidence24/fixtures.json`|`cb706a29d979c0552126baf3107cd8579d55fc95d9af2a94f892b7c2dc2d88e9`|
|`docs/inconsistent-evidence24/srm_cases_INCO24.zip`|`bc569a2ae9db3074d2a124e386a5807072ec955295ed889d8a2833824b40f2bb`|
|ZIP内 `cases/griffiths_1999_coarse/results/griffiths_1999_coarse_critical24.json`|`8d603e36d17aff6fd21dd46adefe40a778364bcbe6c8eab1608860ea835e70ed`|

今回追加した `tests/data/griffiths_2d_reference.json` は、この生結果の節点順・変位と、fixtureの節点ID・座標・拘束をそのまま抽出したもの。推測・補間・文献値からの生成はしていない。元の公開ブックSHAと出典も同梱し、Apache-2.0のライセンスを `tests/data/GRIFFITHS_REFERENCE_LICENSE.txt` に保存した。

JSON内容の固定SHAは `29822140bcfb32da9954d75db23062b65ea28f127630971475ce5c2ab0cb9715`。キーをソートし、空白なし・有限数のみのJSONへ直列化して検査する。Windows checkoutのCRLFは許容するが、変位、点数、節点順、単位、材料、出典の変更を拒否する。任意のSHAを指定して承認するオプションは無い。

従来の文書、既存Rustレポートの手順文字列やF=1評価器にある「661節点データが無い」は、今回の検索前の記述である。最新の所在と手順は本書を参照する。過去の数値・失敗・判定不能はそのまま残す。

## 実測した比較と今回の限界

`tools/griffiths_2d.py archives` は、保存された実求解の全661点について、x/z両成分の非加重L2を奥行き0／0.5／1 mで再計算する。式は `||u3D_xz − u2D||₂ / ||u2D||₂`。元の `displacement_sample_relative_L2 = 0.05` を使い、全奥行きが5％以内のときだけ変位L2を合格とする。最大誤差10％、奥行き差1％も別々に出力する。

|既存ケース|中央面相対L2|5％判定|面内h系列|
|---|---:|---|---|
|xy1-ny1|15.581998707％|不合格|いいえ|
|xy1-ny2|15.885691930％|不合格|いいえ|
|xy1-ny1-global64|15.593219845％|不合格|いいえ|
|material-half1200|15.672391195％|不合格|いいえ|

これは元の661点実データに対する**保存済み結果の再集計実測**であり、新しいFEM求解ではない。奥行き・荷重刻み・材料積分幅の対照をh系列へ混ぜず、連続側は判定不能のまま。再集計JSONは [今回の記録](../interim/2026-11-08-griffiths-2d/README.ja.md) に保存する。

柱のTet4最細重力場6.828％の不合格、連続推定1.034％の合格、Tet10の5％合格は変更していない。新しい斜面h系列の元2D変位5％達成、正確なFs、一般3D空間収束、P6全体を完了とは扱わない。

## ユーザーPCでの準備・求解・比較

元のF=1系列、既存の入力準備・validate・体積220 m³の記録は保持した。既定の `prepare_slope_h_series.py` が生成するmesh／ground／case／parent-quadsの全3段階は起点コミットとバイト一致することを確認した。

公開2D結果は固定低減係数 `1.3500000000000003`、c=`7.407407407407406` kPa、φ=`15.088616807815535`°、E=100000 kPa、ν=0.3、γ=20 kN/m³、ψ=0、非関連流れ、Rankine面なし。**F=1の結果をこの参照変位へ比較することは拒否する。** この低減係数は過去の比較条件であり、物理的な破壊Fsの認定ではない。

明示的な `--reference-2d` を付けると、新しい出力先へこの材料の別系列を準備する。面内1/2/4分割、奥行き1層固定、1200／4800／19200 Tet10、2583／9963／39123節点。全y変位固定、元のx/z拘束、ゼロ種応力から全重力、初期増分1/32・最小1/4096・材料積分幅1e−4を保持する。全661点の前面・中央面・背面を実節点へ直接対応させ、欠落点を補間しない。

以下は**ユーザーPC用**のPowerShell例。このクラウドでは求解コマンドを実行していない。固定依存を用意し、新しい空の出力先を指定する。

```powershell
python -m pip install -r requirements-reference.lock.txt
if ($LASTEXITCODE -ne 0) { throw 'dependency setup failed' }
cargo build --locked --release
if ($LASTEXITCODE -ne 0) { throw 'release build failed' }
$ref = 'tests/data/griffiths_2d_reference.json'
$root = 'runs/griffiths-2d-matched'
python tools/prepare_slope_h_series.py $root --reference-2d $ref
if ($LASTEXITCODE -ne 0) { throw 'preparation failed' }
foreach ($level in @('xy1-ny1', 'xy2-ny1', 'xy4-ny1')) {
    $p = "$root/$level"
    target/release/soil-fem-3d.exe validate "$p/mesh.json"
    if ($LASTEXITCODE -ne 0) { throw 'mesh validation failed' }
    target/release/soil-fem-3d.exe solve-stages "$p/mesh.json" "$p/ground.json" "$p/case.json" "$p/result.json"
    if ($LASTEXITCODE -ne 0) { throw 'solve incomplete; no convergence or Fs verdict' }
}
python tools/griffiths_2d.py series runs/griffiths-2d-comparison.json --reference-2d $ref --folder $root --solver target/release/soil-fem-3d.exe
```

`series` は元の準備器から入力を再生成して、manifestの改変・再ハッシュ、全入力の改変、F=1混入、hや奥行きの変更、欠落点、未完成Run、非有限変位、ソルバーの変更を拒否する。初期釣合いの種からの全変位を比較し、表示リセット後の変位は使わない。力残差・全体力1e−8、正規化降伏1e−7、工学ひずみ0.05、拒否履歴の保持を検査したうえで、既存の `export_stages` と `verify_stage_package` の独立力・材料・生GP/VTU監査を必ず実行する。監査を省略するオプションは無い。入力・結果・ソルバーの実ファイルSHAも記録する。

`five_percent_verdict` は各段階の変位L2のみの判定。`all_levels_five_percent_passed` も離散サンプルの結果である。最大誤差・奥行き差は別ゲートとして必ず確認する。この比較器は場の連続極限を外挿しないため、`continuum_verdict` は `indeterminate`、`griffiths_slope_five_percent_closed` は false を保持する。参照なしRichardsonによる天端・プリズム平均のF=1診断は従来の `assess_slope_h_series.py` に残し、今回の低減済み参照比較へ混ぜない。

終了コードは0＝全離散段階の変位L2が5％以内、1＝比較を実測したが5％不合格、2＝データ欠落・改変・未完成などで比較拒否。拒否時に合格JSONは作らない。既存レポートは上書きしない。監査ログの既存出力も上書きしないので、再評価では新しい結果ディレクトリを用意する。判定不能・不収束・終了コードから破壊側上端やFsを作らない。

保存済み実結果の再集計だけなら、次のコマンドを使う。終了コード1は期待どおりの既存不合格である。

```text
python tools/griffiths_2d.py archives runs/griffiths-archived-2d.json --reference-2d tests/data/griffiths_2d_reference.json
```

## 元アーカイブとの照合の再現

固定版の2ファイルをバイト変更せず取得する。ZIP内の元の変位の節点IDと全値を検査するだけで、Excelや斜面求解は実行しない。

```powershell
New-Item -ItemType Directory -Force work/griffiths-source | Out-Null
$base = 'https://raw.githubusercontent.com/tak063495-prog/excel-vba-soil-fem/440d246d9e32f5892c79539d4a7b797ad29442ee/docs/inconsistent-evidence24'
Invoke-WebRequest "$base/fixtures.json" -OutFile work/griffiths-source/fixtures.json
Invoke-WebRequest "$base/srm_cases_INCO24.zip" -OutFile work/griffiths-source/srm_cases_INCO24.zip
python tools/griffiths_2d.py verify-source runs/griffiths-reference-source.json --reference-2d tests/data/griffiths_2d_reference.json --fixtures work/griffiths-source/fixtures.json --archive work/griffiths-source/srm_cases_INCO24.zip
```

今回の追加試験は、実データ・4不合格の再現、元参照の改変と欠落の拒否、CRLF、既知場の合格／不合格／5％境界、全奥行き・全成分、最大誤差、非有限・ゼロノルム、点数違い、参照条件の準備、直接節点対応、manifestや入力の改変、F=1混入、未完了Run、物理ゲート、監査失敗時に報告しないことを検証する。合成場は単体試験内だけで使用し、実測データとして保存していない。

範囲外は一般3D空間収束、正確なFs、実断面SRM、Native UI、構造・水・JOINT、別PC/GPU、完成配布、乾燥MCの4/8/16・8/16/32と既存MCC掘削274〜14909要素の再求解。これらをP6完了へ読み替えない。
