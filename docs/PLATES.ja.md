# 膜・DKT曲げを持つ平板と地盤・梁の結合

plate-surfaces/3の掘削集合指定を、実Nativeの有限幅二層モデルへ適用した。同じ上層ボディの内部切断面でも保持側を選び、下層の共有ボディ境界と合わせた四周壁72 m²・128板を保持地盤へ結合した。274 Tet10の5工程、掘削途中保存再開、公式GPUの曲げと12 m²の壁断面、生の板値を確認した。従来v1外面・v2共有面の実144要素モデルについてはmesh/case/bindingが旧処理と完全一致し、Python全128件も合格した。一方、945要素の細分モデルは詳細な板属性のNative64 KiB上限と、外部ソルバーの縮約構造結合の記号接続サイズ上限で拒否された。全非局所結合を保つ別作業用修正を検証中で、接触・すべり・一般空間収束の受入と同一視しない。

平板三角形を連続体地盤と梁に接続した。局所軸のCST膜ひずみとDKTの曲げを組み合わせ、節点並進3成分・回転3成分を既存の構造モデルで共有する。構造内部の自由度を静的縮約してMC/MCC共通の全体Newtonへ加える。板は線形弾性、小変形、平面要素であり、厚板の横せん断、材料塑性、座屈、幾何学的非線形を含まない。

DKTの辺接線方向回転をHermiteたわみの辺中央微分から定め、辺法線方向回転は線形にする。膜は平面応力。公開定式化は [Code_Aster R3.07.03](https://code-aster.org/doc/v17/manuals/man_r/r3/r3.07.03/index.html) を参照し、コードは独自に実装した。独立Python参照は頂点・辺条件から12個の二次回転多項式係数を解き、別の5×5 Gauss/Duffy則で剛性を積分する。

面法線回りの回転には、膜のスピンとの差に対する明示的なペナルティを使う。`drilling_factor` は G×板厚に対する係数で、入力必須、許容範囲1e−8〜0.01である。これは曲げ材料則ではない。膜・曲げ・drillingのエネルギーを別々に出力する。この係数と接続部回転の感度を確認せず、特に面内ねじりを受ける梁と板の接続剛性を確定したと扱わない。

## 入力と出力

既存の `beam_system` を共通の構造モデルとして使い、`schema: soil-fem-structures/1`、`plate_sections`、`plates` を指定する。節点IDを梁と共有すると並進・回転を共有する。地盤とは、`match_soil_node` または一節点の `soil_node_set` で並進を結ぶ。地盤連続体の回転付着や任意の面全体への自動結合ではない。

板断面にはID、E、ν、板厚、単位体積重量、drilling係数を指定する。板ごとに3節点ID、断面ID、初期有効状態、一定/比例載荷のグローバル面トラクションを指定する。自重は板厚×面積×単位体積重量から加える。曲面は平板を並べた近似になる。

設置・撤去は `install_structures` / `remove_structures` の `members` に板・梁IDを指定する。既存の `install_beams` / `remove_beams` も同じ処理を呼ぶ。設置時の全並進・回転を部材ごとの基準として保存する。確定/試行状態、切戻し、チェックポイントにも板の有効状態と設置基準を含める。板を含む結果は `soil-fem-stage-result/14` で識別する。

`examples/stages-plate-beam-column.json` と `examples/ground-plate-beam-column.json` は、既存 `hardening-column-tet10-mesh.json` の形状に板2枚・梁1本を接続する弾性対照例である。地盤はMCの弾性枝に留める強度を明示しており、実地盤の強度同定例ではない。板・梁の設置、自重、載除荷、撤去と意図した切戻しを含む。

ルートでrelease build後、次のように計算・出力する。出力先は新しい名前を使う。

```powershell
New-Item -ItemType Directory -Force runs | Out-Null
target/release/soil-fem-3d.exe solve-stages examples/hardening-column-tet10-mesh.json examples/ground-plate-beam-column.json examples/stages-plate-beam-column.json runs/plate-control.json
.venv/Scripts/python.exe tools/export_stages.py examples/hardening-column-tet10-mesh.json runs/plate-control.json runs/plate-control-vtu
```

通常の施工出力に加え、`plates.pvd` は基準三角形と全節点変位・回転、`plate-samples.pvd` / CSVは三つの生の曲率・膜力・曲げモーメントを出力する。局所軸と単位を保持し、節点への応力平滑化を行わない。非有効時は空フレームを明示する。SolidCraft内の板結果専用レイヤーは未接続である。

## 受入の範囲

Rustの剛体6モード・一定膜ひずみ/曲率・エネルギー・回転座標変換・面荷重の力/モーメント、薄板片持ちの理論解を検査している。MC/MCC施工への接続、保存再開、独立全体3D計算、VTU/生値出力の最終結果は固定ログで確認する。試験開始や部分合格だけでP3b/P6全体を完了にしない。

`runs/plate-foundation-20261026-v1/control/verification.json` は548ソースと実行版 `8a3eda114a0d92a46f5854b84e08c010928c513f7b2a0ccf0fc74dfd86da2e7d` を固定した限定検証である。12 Tet10・324地盤積分点・板2枚・梁1本、5フレーム・24受理増分・2回の意図した切戻しを完了した。全節点地盤変位の独立相対差は最大2.456e−12、全GP応力は2.767e−12。構造節点の全並進/回転、設置基準、板の全評価点の曲率・膜力・曲げモーメント・エネルギー、梁端力も比較した。参照解は元入力から独立に組み立て、Rustの計算済み変位・応力・設置基準を生成に使わない。

初期化、設置、載荷、除荷、撤去の5地点から別CLIへ保存再開し、最終状態・全フレーム・施工監査・時間以外の試行履歴が連続計算と一致した。意図した切戻しを除いた経路でも最終状態と全フレームが一致した。地盤の塑性履歴/塑性仕事は全点で0の対照である。板と生評価点のVTUも再読込して全フィールドが一致した。これは塑性地盤の一般壁変位や空間収束の受入ではない。

別のMC/MCC柱に板を設置・載除荷・撤去する試験を含め、専用Rust試験4件はrelease版で合格した（`runs/plate-foundation-20261026-v2-focused.log`）。総合Rustも247件・無視0で完了し、`runs/plate-foundation-20261026-v1/full-rust-terminal.json` にソース/実行版を結び付けた。初期のPython55件と全ターゲットclippyに加え、後続の面対応3試験を含む全Python58件も合格した。初回のdebug試験ではMCC初期荷重用のtop面が試験メッシュに欠けていたため拒否された。面セットを補い、ソルバーの入力検査を変更せず再検証した。

後続でSolidCraft外表面から平板へ対応させる処理と、粗密2メッシュ・drilling3条件・独立全場照合を追加した。[CAD壁面の対照例](CAD_PLATE_WALL.ja.md)を参照。共有節点間の補間の違いを残件として明示し、山留め壁や一般収束の受入へ広げない。

板の一般的な空間収束、drilling依存性、土留め壁・切梁の実CADでの全体比較、地盤界面の分布伝達とJOINT、任意面からの入力、全水条件との組合せ、専用の画面表示は残る。P1体積投影・3場混合MCCの部材接続も未対応であり、既存の拒否を保持する。全体の未達項目は [P6一覧](P6_STATUS.ja.md) に保持する。


## 掘削集合で保持側を選ぶ同一ボディ内の壁（開発版 /3）

`tools/plate_surfaces.py` の `soil-fem-plate-surfaces/3` では、同じ材料ボディ内部を切る掘削集合に壁を対応させられる。`attachment: "shared_excavation_interface"` と `opposite_excavation_set` を明示する。面ごとに、指定した `owner_body` に属して掘削集合の外にある保持要素が一つ、集合内の掘削要素が一つ隣接することを要求する。ボディ名が両側で同じ場合も、この集合から向きを確定する。

```json
{
  "schema": "soil-fem-plate-surfaces/3",
  "sections": [{"id":"wall","young_kpa":30000000,"poisson":0.2,"thickness_m":0.2,"unit_weight_kn_m3":24,"drilling_factor":0.0001}],
  "surfaces": [{
    "id":"pit_wall", "face_set":"Pit:wall", "owner_body":"Soil",
    "attachment":"shared_excavation_interface", "opposite_excavation_set":"Pit",
    "section_id":"wall", "initially_active":false
  }]
}
```

`face_set` は入力メッシュ上に用意する。底面を含めず壁面だけを選ぶ場合も、呼出側で面集合を分ける。`opposite_body` との重複指定、空・不正・重複した要素集合、保持側がない面、集合の内部同士の面を拒否する。Tri6の全中間節点、両側の同一節点と逆向き法線を検査し、保持側の向きで四つのTri3板へ分ける。元の座標・接続・入力を変更せず、生成した構造節点と部材群を返す。

既存の /1 の外面、/2 の二ボディ共有面の出力契約は保持する。接触・すべりや、CAD編集後の永続的な面対応を追加する機能ではない。保持側と掘削側の共通並進、板の回転自由度、二次地盤と区分線形板が節点間では同一の補間にならない制約も従来どおりである。同一ボディの保持側支持、掘削集合の反転に伴う向き、集合不備、掘削側中間節点不整合を4試験で検査し、既存7件を含む11件が合格した。有限幅の実CAD・力学・独立解との受入は継続する。
