# 地盤に結合する3D梁・切梁

2026-10-24。小ひずみの線形3DフレームをMCC/MC全体FEMへ接続した。軸伸び、2方向曲げ、ねじり、Euler–Bernoulli/Timoshenko断面、任意の局所軸、回転端部解放、分布自重、節点力/モーメント、初期軸力、設置/撤去を扱う。梁の塑性・座屈・大回転、引張専用アンカーは[別モジュール](ANCHORS.ja.md)へ接続済み。板/シェル、履歴JOINTは別要件であり、P3b全体の完了とはしない。

## 操作

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-Stages.ps1 -Design examples\geostatic-slope.json -GroundFile examples\ground-beam-slope.json -StageFile examples\stages-beam-slope.json -ElementOrder 2 -MeshSizeM 4
# 塑性MCC地盤と同じフレームの結合
.\scripts\Invoke-Stages.ps1 -Design examples\geostatic-slope.json -GroundFile examples\ground-geostatic-mcc-slope.json -StageFile examples\stages-beam-slope.json -ElementOrder 2 -MeshSizeM 4
# 理論・独立連立解・切戻し・保存再開・既存掘削の回帰を確認
.\scripts\Test-Beams.ps1
```

地盤はSolidCraftの閉ソリッドから取得する。梁の節点・中心線・断面・荷重・設置工程は現段階では解析ケースの`beam_system`へ入力する。SolidCraft画面内の梁編集パネルや中心線の自動読込はP5の残作業。地盤と梁は同じm–kN–kPa座標系を使う。

## 入力と接続

`beam_system.schema`は`soil-fem-beams/1`。

|項目|意味|
|---|---|
|`nodes[].id`, `position_m`|梁節点名とCAD基準位置。梁だけの自由節点も置ける|
|`soil_node_set`|ちょうど1節点を含む地盤節点セットへ並進3成分を接続。座標一致を検査|
|`match_soil_node: true`|座標が一致する地盤節点へ接続。候補0個/複数個を拒否。`soil_node_set`と併用しない|
|どちらも未指定|自由な梁節点。同じ場所に地盤節点があっても自動では結合しない|
|`sections`|E[kPa]、ν、A[m²]、Iy/Iz/J[m⁴]、単位体積重量[kN/m³]|
|`formulation.kind`|`euler_bernoulli`または`timoshenko`。後者は有効せん断面積Ay/Az[m²]を指定|
|`members`|節点i/j、断面ID、局所y方向、初期有効状態|
|`rotational_releases[2][3]`|両端の局所rx/ry/rz解放。等価分布荷重も同じ解放条件で処理|
|`initial_axial_force_kn`|引張正の初期軸力。梁の初期応力であり、片効きアンカーのモデルではない|
|`constant_line_load_kn_m`, `boundary_line_load_kn_m`|全体座標の一定/載荷係数比例分布荷重。断面重量×Aの鉛直自重も加える|
|`nodal_loads`|全体座標の節点力[kN]とモーメント[kN·m]。一定分と`boundary`係数比例分|
|`translation_constraints`, `rotation_constraints`|方向ベクトルに沿う値。並進[m]、回転[rad]。地盤に接続した並進は地盤側境界で拘束する|

梁回転は独立したDOFである。地盤への接続は共有並進のみで、地盤の回転勾配との完全付着ではない。明示した回転拘束の反力モーメントは外部支持へ集計し、地盤へ伝えたものと表示しない。地盤要素内部への補間接続や回転の面結合は未実装。削除した地盤節点に有効梁を接続したままにすると入力エラーになり、先に梁を撤去する必要がある。

## 求解と釣合い

部材ごとの軸・曲げ・ねじり剛性を全体座標へ変換し、解放端の回転を厳密に消去する。自由梁節点の並進と回転を含む内部系をCholeskyで解き、地盤との接続並進へSchur補行列を戻す。地盤側CSRには離れた接続節点間の非局所成分も組み込む。

`rotation_length_scale_m`で回転未知量の単位を調整し、角度と並進を同じ数値尺度で扱う。`maximum_condensed_dofs`は内部密行列の上限で、最大3000。大規模フレームの疎ブロック求解は後続作業。剛体ねじり・拘束不足・解放で生じた未使用回転を擬似ばねや対角シフトで固定せず、入力エラーにする。ベンチマークでは解放により未使用となる端部回転を明示拘束し、その反力がゼロになることを確認する。

内部系の代数残差に加え、回復した物理節点力/モーメントを自由方向へ投影した真残差を検査する。梁の外力、地盤接続へ含めた外力、外部支持反力を区別し、地盤と梁を合わせた力・原点まわりのモーメントを集計する。地盤の接続節点力に含まれる梁外力を二重計上しない。

## 工程と状態

```json
{"kind":"install_beams","name":"install_frame","members":["Left","Right"],"placement":"current_displacement"}
{"kind":"boundary","name":"load_frame","target_factor":1.0}
{"kind":"boundary","name":"unload_frame","target_factor":0.0}
{"kind":"remove_beams","name":"remove_frame","members":["Left","Right"]}
```

設置時の12成分の並進・回転を各部材の参照状態へ記録する。配置は`current_displacement`または`cad_new_nodes`を明示する。後者は新たに有効となる梁自由節点の保存変位/回転をCAD基準の0へ戻す。地盤に接続した並進は常にその時点の地盤変位を引き継ぐ。既存の有効梁節点は置き直さない。

部材の自重と初期軸力を加えた設置直後の内力を一旦保持し、地盤への保持荷重を段階的に解放する。撤去も同じように荷重再配分を段階的に進める。梁内部の自由DOFは各試行で平衡を解く。部材長・局所軸・剛性はCADの弦に基づく固定値で、設置後形状を使う幾何学的非線形ではない。

有効部材、設置参照12成分、全梁節点の物理変位/回転を確定状態に保存する。切戻し・再開では地盤の全積分点履歴とともに復元する。チェックサムを再計算したデータでも、回復した梁変位/回転が内部平衡と一致しなければ再開を拒否する。撤去後の値は非有効の保存値で、部材力は生じない。

梁がある場合の結果スキーマは`soil-fem-stage-result/7`。梁なしの既存結果の数値と省略フィールドは維持する。

## 出力

- `beam_results`：節点変位/回転、局所軸、両端作用力、解放端の相対回転、始端/中央/終端の断面力、支持反力、真残差。
- `beams.pvd`と`beams-*.vtu`：有効梁だけの線形状。初期未設置・全撤去の時刻も空のVTKを明示し、前時刻の梁を残さない。
- `beam-members.csv`：N/Vy/Vz[kN]とT/My/Mz[kN·m]。Nは引張正。節点端部作用力と断面力の符号を区別する。
- `beam-nodes.csv`：CAD位置、全物理変位、全回転、有効状態。`beam-fields.json`にIDと表示の意味を保存する。

梁の`total_displacement_from_cad_m`は地盤の表示RESETと独立した物理変位。地盤の施工表示変位と混ぜて同じ基準の値としない。VTU形状はCADの直線弦なので、この全変位をWarp By Vectorに使う。梁の内部曲げ曲線を細分して表示する機能は未実装。

定式化の参考は[TU DelftのEuler–Bernoulli梁](https://teachbooks.tudelft.nl/computational-modelling/structural_linear/euler_bernouilli.html)、[Timoshenko梁](https://teachbooks.tudelft.nl/computational-modelling/structural_linear/timoshenko.html)、[OpenSeesの3D梁・端部解放の仕様](https://openseesdocumentation.readthedocs.io/en/latest/user/manual/model/elements/elasticBeamColumn.html)。ソースの複写ではなく独立実装。Timoshenkoは静的な断面柔性と整合する閉形式で、単純な等次線形補間の完全積分によるせん断ロッキングを持ち込まない。

実測と適用範囲は[梁の検証記録](BEAMS_VERIFICATION.ja.md)、目標全体の残りは[P6受入一覧](P6_STATUS.ja.md)で管理する。
