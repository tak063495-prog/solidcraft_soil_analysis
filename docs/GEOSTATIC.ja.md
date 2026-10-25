# 一般地形の初期応力・重力平衡と施工変位

2026-10-24。水平直方体多層の`horizontal_k0`に加え、`equilibrate`経路を実装した。SolidCraftのボディ名と材料を対応付け、初期有効応力と材料履歴を読み込み、自重・既知水圧・外力・拘束に対する非線形平衡を求める。斜面へK0を代入しただけの状態を施工の起点にはしない。

初期化と施工は同じTet4/Tet10、全積分点MCC/MC、Newton・ラインサーチ・疎求解・切戻しを使う。初期化中に形成した塑性ひずみ、MCCのpc・間隙比、MCの履歴を施工へ引き継ぐ。有限変形、地形の自動復元、地質年代の堆積史の同定は含まない。

## 実行

```powershell
cd D:\Work\soil-fem-3d
# SolidCraftの実斜面形状を取得し、初期平衡→載荷→除荷
.\scripts\Invoke-Stages.ps1 -Design examples\geostatic-slope.json -GroundFile examples\ground-geostatic-mcc-slope.json -StageFile examples\stages-geostatic-slope.json -ElementOrder 2 -MeshSizeM 4
# 積分点ごとのK0種応力→無載荷斜面の平衡→CAD指定の局所掘削
.\scripts\Invoke-GeostaticSlopePit.ps1
# 理論/独立計算、保存再開、既存MCC結果の回帰を含む一式
.\scripts\Test-Geostatic.ps1 -MccBaseline runs\tet10-boundary-test-20261023-224400-66dd339e\flat\stage-results.json
```

最初の例の100kPaの初期面圧は明示した検証用拘束荷重である。2番目の例は初期地表面圧0で、単一の平面状地表`z=10−x/3`の上載土被りから生成する例。既定は種状態K0x=K0y=1、pc/p=4の過圧密ベンチマークで、自動推定値ではない。`-K0`と`-PcOverMean`で明示的に変更できる。水平床z=7m、平面4×4m、平均深さ1mの掘削を扱う。現場一般の地圧を推定する機能ではない。

初期平衡だけを作る場合は`Invoke-Ground.ps1`へ`-InitializationCaseFile`を指定する。拘束、初期境界係数、Newton設定をそのケースから読み、施工工程を空にして平衡化する。平衡・独立組立・出力検査に合格した結果だけに`accepted-run.json`を作る。

## 入力の意味

`initialization.method: equilibrate`では`layers`の代わりに`bodies: [{body, material_id}]`を指定する。MCCの材料は従来どおり`mcc`、MCは`mohr_coulomb`を使う。一般形状のボディには水平層の上下標高や架空のOCRを要求しない。

`initialization.geostatic`の主な項目：

|項目|意味|
|---|---|
|`seed.kind: body_states`|ボディ名をキーにした均一な初期材料状態|
|`seed.kind: integration_points`|メッシュに結び付いた全積分点の初期状態|
|`initial_boundary_factor`|初期平衡完了時の一般拘束値・追加荷重の係数|
|`maximum_initial_engineering_strain`|初期化の種状態から許す工学ひずみ成分の絶対最大値。例は0.05|
|`relative_volume_tolerance`|積分したボディ体積と入力CAD体積の許容相対差|
|`base_loads`|施工中も残す節点荷重・面圧・トラクション|
|`hydrostatic_face_sets`|初期外表面のうち外水圧を与える面セット。新しい掘削面は別途自動処理|
|`verification_reject_once`|検証専用の収束試行棄却。通常の解析では省略|

応力は有効・圧縮正、並びは`xx,yy,zz,xy,yz,zx`。塑性ひずみも圧縮正。履歴の型、有限値、材料降伏条件を検査する。MCCは正の平均有効応力とpc、間隙比とNCL/除荷線の整合を要求する。MCへpcを付けない。種状態で降伏面の外にある応力を自動で切り詰めて受理する処理はない。

積分点入力では`mesh-identity MESH.json NEW_IDENTITY.json`で出力する正規化メッシュSHA256とレイアウトを使う。各点へ`element_index, local_index, material_id, position_m, state`を保存する。点数、順番、位置、材料、メッシュSHAの不一致を拒否する。要素次数や積分次数を変更したら再生成する。

形状を先にメッシュ化する場合は`Invoke-Ground.ps1 -PrepareOnly`を使う。これは`prepared-model.json`に`equilibrium_verified: false`を記録する。積分点入力作成後、`Invoke-Stages.ps1 -PreparedModel <folder> -GroundFile <new physics> -StageFile <case>`で明示した新しい物性入力へ置き換える。元のSTEP/メッシュ/物性のハッシュと置換を記録し、平衡済みパッケージには偽装しない。

## 平衡化と履歴

種状態の節点内力を`I_seed`、種状態の実荷重を`F_physical(0)`として、仮保持荷重を`H=I_seed−F_physical(0)`とする。進行度λで`F_physical(λ)+(1−λ)H`を解き、λ=1で仮保持荷重をゼロにする。既存の不平衡を一度に除く代わりに、実際の材料更新を伴う荷重経路として段階的に解放する。最終的に物理荷重だけで自由度残差を再検査する。

これは指定した種状態からの平衡経路であり、初期応力が一意に定まる保証ではない。地質履歴・測定応力・境界を選ぶ責任は入力側にある。種状態、仮保持荷重、全試行、受理状態、塑性変化、初期ひずみ、実荷重残差を保存し、どの経路から形成したかを追跡できる。

初期ひずみ上限を超えた収束試行は確定しない。未完了結果に最終確定履歴と失敗理由を残し、施工へ進めず、受入済みVTUを作らない。形状を更新して変形量を隠す操作も行わない。

## 施工・掘削と出力

初期平衡時の変位を`reference_displacement_m`に保存する。内部の累積変位・応力・全材料履歴はそのまま保持する。監視値とVTUの`displacement_m`は「累積変位−初期平衡変位」で、施工に伴う沈下/隆起を示す。初期化変位と種状態からの累積変位は別フィールドで出力する。一般初期化の結果スキーマは`soil-fem-stage-result/4`。

`maximum_absolute_engineering_strain_from_seed`と`...from_initial_equilibrium`も区別する。初期平衡化が大きく変形したかと施工増分の大きさを別々に調べる。メッシュは小ひずみの基準形状を維持する。

局所掘削はSolidCraftの閉じたツールボディとの実Boolean交差で要素を選ぶ。地表の傾斜面から掘り始めることができる。周辺地表セットには残存する上向きの元外表面を使う。底面支持から切り離された残存地盤や、露出面につながらない埋没空洞は拒否する。現在の支持連結検査は水平なモデル底面を前提とする。

掘削前の確定応力から、残存地盤に働く保持荷重を計算し、重力・外力の撤去と切断面の応力解放を段階的に進める。初期化の仮保持荷重と掘削の解放荷重は別々に保存する。元表面の`base_loads`は撤去面では消え、新しい切断面へ自動移動しない。撤去される節点への固定節点荷重は入力エラーになる。

地形の面圧は`free_surface`等の一般面荷重で指定する。従来のスカラー`surface_surcharge_kpa`/pressure工程は水平面の機能を維持しており、斜面全体の上載荷重を自動的に意味しない。`free_surface`は外向き法線のz成分が正の面であり、CAD面の永続IDや任意の現場地表分類を代用するものではない。

チェックポイントは初期化途中、初期化完了直後、施工途中を扱う。全履歴と初期化/掘削保持荷重、表示基準、増分経路、入力/バイナリSHAを復元する。別工程ケースへ`GroundRun`を渡すと種状態からそのケースで再平衡化する。既存計算の形成履歴をそのまま継続する場合は同一入力の`resume-stages`を使う。

## 現在の範囲

固定水位の静水圧と、明示した初期外水面・新しい掘削面の水圧に対応する。水位が境界三角形を横断する場合はCAD分割を要求する。水位変更、部分水没積分、Darcy/Richards、自由水面、圧密時間は未実装。非線形MCCの一般3D塑性変位を独立した別ソルバーで全面照合した段階ではない。検証の内訳は[検証記録](GEOSTATIC_VERIFICATION.ja.md)を参照。P6全体の残りは[P6一覧](P6_STATUS.ja.md)で管理する。
