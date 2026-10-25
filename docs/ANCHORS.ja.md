# 引張専用アンカーと埋込み接続

2026-10-24。直線状の弾性・引張専用アンカーをMCC/MC全体FEMへ接続した。自由長、分割した付着長、力制御による緊張、平衡後のロックオフ、弛緩/再緊張、設置/撤去を扱う。付着部は**完全付着**で、引抜き耐力や付着すべりを表すモデルではない。

## 実行

```powershell
cd D:\Work\soil-fem-3d
# MC弾性地盤との独立全体比較用
.\scripts\Invoke-Stages.ps1 -Design examples\geostatic-slope.json -GroundFile examples\ground-beam-slope.json -StageFile examples\stages-anchor-slope.json -ElementOrder 2 -MeshSizeM 4
# 塑性MCC地盤での同じアンカー工程
.\scripts\Invoke-Stages.ps1 -Design examples\geostatic-slope.json -GroundFile examples\ground-geostatic-mcc-slope.json -StageFile examples\stages-anchor-slope.json -ElementOrder 2 -MeshSizeM 4
.\scripts\Test-Anchors.ps1
```

地盤形状はSolidCraftから取得する。アンカーの直線・断面・材料・工程は現段階では解析ケースに入力する。画面内の中心線作成/材料/工程パネルはP5の残作業。

## 自由長・付着長と地盤への接続

`anchor_system.schema`は`soil-fem-anchors/1`。各アンカーに`id`、E[kPa]、A[m²]、`head`、`bonded_path`、`initially_active`を与える。`head`から付着経路の最初の点までが自由長。付着経路の各区間が分割した付着長となる。点は直線上で順番に並べる。自由長・各付着区間長はCAD座標から計算し、区間の軸剛性はEA/Lとする。

```json
{
  "schema":"soil-fem-anchors/1",
  "embedding_tolerance_m":1e-10,
  "anchors":[{
    "id":"Tie1", "young_kpa":200000000.0, "area_m2":0.0005,
    "head":{"position_m":[0,4,9],"body":"SoilSlope"},
    "bonded_path":[
      {"position_m":[4,4,6],"body":"SoilSlope"},
      {"position_m":[5.5,4,4.875],"body":"SoilSlope"},
      {"position_m":[7,4,3.75],"body":"SoilSlope"},
      {"position_m":[8,4,3],"body":"SoilSlope"}
    ],
    "initially_active":false
  }]
}
```

上の例は自由長5m・付着長5m。自由区間の途中を地盤へ付着させず、頭部と付着開始点だけを結ぶ。付着区間は各点の変位を地盤の形状関数から補間し、隣り合う区間張力の差をその点で地盤へ渡す。頭部は地盤または連続体の壁に埋め込める。自由な梁節点への直接接続は未実装。

各点はTet4/Tet10ホストを座標から探索する。曲面Tet10は逆写像を解き、Bernstein制御点の境界箱、自然座標の領域内判定、座標再現、分割の1を検査する。共有面では同じ補間となる全ホストを保持し、少なくとも1つが残れば接続を維持する。非整合/重複ホストが違う変位補間を与える場合は拒否する。`body`を与えると指定ボディの中だけを探索でき、亀裂/接触面の反対側を誤って選ぶことを防ぐ。

力は補間の転置でホスト節点へ渡し、力・モーメント・仮想仕事を保存する。Tet10の負の形状関数をゼロへ切り詰めない。最近傍節点1点への荷重集中や、未登録の外部点をモデル内へ投影する処理はしない。有効アンカーのホストを掘削で全て撤去する入力は拒否し、先にアンカーを撤去する。

地盤の体積から鋼材体積を差し引かない追加剛性の埋込みモデルである。アンカーは質量なしで、曲げ・自重・熱・クリープ・鋼材塑性・破断は扱わない。付着すべり/引抜きを評価する場合は伝達則を追加する必要がある。

## 引張・弛緩と緊張工程

固定CAD軸tに沿う区間伸びを`d = t·(u_j−u_i)`とし、定着後の張力は`N = max(0, EA/L × (d−d_ref))`。圧縮力は生じない。張った側はEA/L、弛緩側は0の接線で全体Newton/ラインサーチを解く。軸や長さを変形後の値へ更新する大変形モデルではない。

伸びが参照と等しい無応力境界では、張った側の接線を使う。`segment_taut`はこの境界を含み、張力0でも参照より短ければ弛緩側になる。

初期張力の幾何剛性（P–Δ）は含めない。アンカーだけで底面から切り離した土塊を支持する構成は、現在の地盤支持連結検査の対象外。

```json
{"kind":"install_anchors","name":"install_tie","anchors":["Tie1"]}
{"kind":"jack_anchors","name":"jack_tie","anchors":["Tie1"],"target_force_kn":50.0}
{"kind":"boundary","name":"preload","target_factor":1.0}
{"kind":"boundary","name":"unload","target_factor":0.0}
{"kind":"jack_anchors","name":"detension_tie","anchors":["Tie1"],"target_force_kn":0.0}
{"kind":"remove_anchors","name":"remove_tie","anchors":["Tie1"]}
```

設置はその時点の全区間伸びを参照へ保存し、無応力で行う。緊張工程は自由区間を現在張力から指定張力へ力制御で移行させる。付着区間の剛性と地盤の平衡は同時に解く。最終的に自由区間が指定張力Pで平衡した変位dから`d_ref = d−P L/EA`を記録し、力制御を解除する。以後は地盤変形に従って増減する張力になり、Pを外力として重ねて加えない。`last_lock_off_force_kn`は過去の記録で、現在張力を固定する値ではない。

緊張解除も力制御0まで進めて定着する。付着部に残る張力は別に計算し、自由区間の張力0だけを理由に付着部の履歴を消さない。撤去はその時点のアンカー内力を保持荷重で徐々に解放する。

既に弛緩している自由区間を0kNで再定着するときは、現在の伸びと既存参照の大きい方を使い、残っている弛みを保つ。正の緊張力を与えずに自然長を短くしない。

自由長の力制御、最終平衡、制御解除という工程の参考は[Itascaのケーブル要素説明](https://docs.itascacg.com/flac3d700/common/sel/doc/manual/sel_manual/cables/cables.html)。同製品の付着すべり則や実装コードは使用していない。本実装の完全付着と、すべり抵抗を持つグラウトモデルを区別する。埋込み補間・仮想仕事の考え方は[embedded element法の研究](https://link.springer.com/article/10.1186/s40323-022-00223-x)も参照する。

## 保存と出力

結果スキーマは`soil-fem-stage-result/8`。設置参照、現在の力制御値、過去のロックオフ力、全区間の張力/張弛状態を確定履歴へ保存する。緊張途中のチェックポイントでは力制御を保ち、工程終端では解除して再開する。切戻しで試行の張力や参照が漏れないことを検査する。チェックサムを再計算した張力キャッシュ改変も物理再計算で拒否する。

- `anchor_results`：点の位置/全変位/地盤へ渡す力、区間長・伸び・張力・張弛状態、弾性エネルギー、ホスト/重み。
- `anchors.pvd`、`anchors-*.vtu`：自由区間と付着区間、有効アンカーだけの線形状、張力と地盤伝達力。未設置/全撤去時刻も空のVTKを明示する。
- `anchor-segments.csv`：全区間の伸び・参照・張力・張弛・自由区間/付着区間の区別。

形状はCADの直線、Warp用の`total_displacement_from_cad_m`は物理全変位で、地盤表示RESETとは別基準。点の地盤伝達力はアンカーから地盤へ作用する向きで、FEM内部力とは逆符号である。補間による節点力とアンカー初期張力を二重に載荷しない。

操作結果と数値は[アンカー検証記録](ANCHORS_VERIFICATION.ja.md)、残る全体要件は[P6一覧](P6_STATUS.ja.md)を参照。
