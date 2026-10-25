# 曲面Tet10と一般境界条件

SolidCraftの地層と掘削ツールを共形分割し、10節点四面体の3変位成分を弾性解析とMCC/MC全体Newton解析へ渡す。水平直方体層のK0経路と、一般形状の初期応力読込・重力平衡経路を持つ。後者の追加仕様は[一般初期平衡](GEOSTATIC.ja.md)を参照。

## 実行

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-LocalExcavation.ps1 -MeshSizeM 4
.\scripts\Invoke-LocalExcavation.ps1 -Design examples\local-excavation-cylinder.json -MeshSizeM 4
.\scripts\Invoke-LocalExcavation.ps1 -QuadratureOrder 4 -MeshSizeM 4
.\scripts\Test-Tet10AndBoundary.ps1
```

`Invoke-LocalExcavation`の既定は`-ElementOrder 2 -QuadratureOrder 3`。`Invoke-Ground`と`Invoke-Stages`で新しくCADから作るときの既定は従来どおり1次。`-GroundRun`を指定した施工解析は保存済みメッシュの次数を使用する。Tet4比較は`-ElementOrder 1`。次数の増加は自由度・積分点数・Newton計算時間を増やす。

## メッシュと積分点

`soil-fem-mesh/3`は既存の四隅`tetrahedra`に加え、`tetrahedra10`、`quadrature_order`を持つ。中間節点は辺01,12,20,03,13,23の順。Gmshの参照節点座標からこの順に並べ替え、VTKのTet10順に統一する。Gmshで曲面へ投影された中間節点を保持し、直線辺の平均位置へ置換しない。

3次元Duffy変換とGauss–Jacobi積から正重み27/64/125点を生成する。27点則は参照座標の次数5以下の多項式を正確に積分する。曲面上の剛性や非線形材料の積分は一般に多項式ではないため、積分則の感度確認も必要。

Tet10のdetJは参照座標の3次多項式。Bernstein係数の正値という十分条件を調べ、判定できない場合は参照四面体を再分割する。負値、ほぼ特異、深さ/計算予算内に正値確認できない要素は拒否する。これは浮動小数点による保守的品質検査で、区間演算による数学的保証とは区別する。積分点だけで見逃す反転を試験に含める。共有辺の中間節点不一致、孤立節点、頂点との取り違えも拒否する。

全体の状態配列は要素順・局所積分点順で格納し、位置、重み、形状関数を`integration_points`へ保存する。MCCの応力、pc、間隙比、塑性ひずみ、試行分割数を全積分点で保持する。撤去要素の全点履歴を凍結し、切戻し・チェックポイント再開でも配列全体を復元する。

Tri6面荷重は曲面の面積ベクトルと6個の形状関数で整合積分する。自重もNiに対して体積積分する。掘削時の保持荷重`H=F_old−F_new−F_int,removed`と段階解放式はTet4と共通で、曲面の切断面にも適用する。

SolidCraftの表示体積は細分三角形による近似値のため、円柱ではSTEPのB-rep体積とわずかに異なる。解析用`volume_m3`はOpenCASCADEのSTEP体積とし、`reported_tessellated_volume_m3`と差率を別保存する。両者の差が0.2%を超える入力は拒否する。Boolean前後の地層体積の検査を緩める意味ではない。曲面のメッシュ体積/面積は引き続きCADからの離散化誤差を持つ。

## 一般境界の入力

`mesh.node_sets`は節点番号集合、`mesh.face_sets`は三角形の四隅側3節点で識別する面集合。Tet10では所有要素からTri6の向きと中間節点を復元する。自動生成集合は`xmin/xmax/ymin/ymax/bottom/top`と`Pit:floor/wall/surrounding_surface`。CAD編集後も同一の面を追跡する永続面IDは未実装。

`examples/mechanics-pressure.json`は弾性用入力。次のコマンドで一般弾性を実行できる。

```powershell
target\release\soil-fem-3d.exe solve-mechanics MESH.json examples\mechanics-pressure.json NEW_RESULT.json
.venv\Scripts\python.exe tools\export_mechanics.py MESH.json NEW_RESULT.json NEW_OUTPUT_DIRECTORY
```

`boundary.constraints`は`direction · u = value_m + gradient · (x−reference_point_m)`を節点集合へ指定する。directionは正規化し、単位方向に換算した値で整合性を検査する。各節点の直交零空間Tとオフセットu0を作り、`u=u0+Tq`として`TᵀKT`を組み立てる。矛盾する拘束、未拘束の剛体モードを持つ連結成分を拒否する。回転自由度や節点間の一般MPCを表す入力ではない。

`nodal_loads.total_force_kn`は集合全体の合力であり、節点数で均等配分する。`surface_loads`は`kind: pressure`または`traction`で、圧力は正が内向き、tractionは全体座標のkPaベクトル。位置による一次勾配も指定できる。面積、曲面法線、Niを用いて節点へ積分する。有効地盤の外面でない面への有効荷重は拒否する。

## MCC施工への接続

`examples/stages-prescribed.json`は`boundary`と`kind: boundary, target_factor`を使い、拘束値と追加荷重を共通係数で増減する。`side_rollers`では既存ローラーへ一般拘束を追加し、`explicit`では一般拘束だけを使う。`oedometer`と一般境界の併用は拒否する。K0初期状態が新しい拘束系で平衡することも必要。

係数0でも拘束方向は保持し、値が0になる。拘束そのものの解放・方向変更を表す操作ではない。強制変位を0へ戻す工程では残留反力が生じ得るため、荷重除去後の自由な残留沈下と区別する。初期MCC状態からの変位、全荷重・拘束係数を確定状態に保存し、試行開始時から拘束オフセットを満たす変位で材料更新する。拒否時はオフセット・積分点履歴も復元する。

## 結果と適用範囲

`results/stages.pvd`は有効地盤の変位・体積平均状態、`results/integration-points.pvd`は生の積分点状態。`volume_average_`の場は重み付き要素平均であり、局所塑性域の最大値ではない。`stages.csv`の監視平均は節点算術平均。正のuZは隆起、負は沈下。周辺に両方が出る場合は最小/最大も確認する。

異なる積分則の独立弾性計算、MCC弾性枝の独立3D速度式、塑性MCCの独立1D基準、全履歴の切戻し・別プロセス再開を組み合わせて検証する。局所塑性3D解全体の独立解、ほぼ非圧縮ロッキングの解消、空間収束、実地盤の精度を証明するものではない。一般斜面の初期平衡、支保・接触、水位変化/排水/浸透、有限変形は後続課題。
