# 将来の本体出力との比較契約案

本体の出力機能の追加ではなく、独立スクリプトが読む交換形式の案。B-03の受入仕様を決定しない。付属 `results/synthetic-reference.json` / `.csv` は人工の90°折れ目で、本体出力ではない。

`geometry_reference.py compare trusted-reference.json body-output.csv report.json` の第1引数には、入力の元の面・節点・向き付き法線・e・tから独立に作った固定参照を与える。比較対象の出力から参照値を作らない。人工fixtureを実モデルへ流用しない。

JSONは `metadata` と `rows` だけを持つ。metadataは次のとおり（追加・欠落も拒否）：

```json
{
  "schema": "b03-geometry-v1",
  "length_unit": "m",
  "coordinate_basis": "global_cartesian_right_handed",
  "normal_sign": "oriented_face_winding",
  "offset_sign": "attachment=neutral+e*normal"
}
```

各rowは `shared_node_id`、`face_id`（空白なしの非空文字列）、`neutral_m`、`normal`、`e_m`、`t_m`、`attachment_m` のみ。3つのベクトルは3成分。欠損ID、重複する節点/面の組、参照にないID、余分な/欠損フィールド、非有限値、非正厚さ、非単位法線、範囲外eを拒否する。一つの面IDでは法線・e・tを一定とする、区分的平面の契約である。曲面の可変法線への拡張は含まない。行の並び順は自由で、IDの組で対応付ける。

`neutral_m` は既存文書のx_p、`attachment_m` はx_s。既存の `position_m` は付着点なので、これを中立点へ読み替えない。面の節点順は別の固定入力記録に残し、その外積から参照法線を定める。法線だけから本来の外向き/内向きは決まらないため、信頼できる向きの記録がない場合は「未確認」とし比較を進めない。

CSVのヘッダーはfixtureの1行目と同じ順序。metadata5列を各行に繰り返し、ID2列、e_m、t_m、neutral_m_x/y/z、normal_x/y/z、attachment_m_x/y/zの順とする。JSON/CSVともUTF-8、数値は丸めず倍精度の往復表記を使う。スクリプトは単位変換、座標基底回転、法線の自動正規化/反転、eの絶対値化・クリップ、付着点の平均化をしない。

1. metadataの完全一致、全ID、有限性・厚さ・法線ノルム・範囲を検査する。
2. 各成分を元入力から作った固定参照と厳密数値比較する。法線とeを同時反転し、e nだけを一致させても拒否する。クリップ後のeが範囲内でも元参照と違えば拒否する。参照自体が範囲外eならその参照も拒否する。
3. 付着点の恒等式 `attachment=neutral+e*normal` の残差を全行保存する。丸め差を含む非ゼロ差は許容せず、契約不一致として拒否する。

新しい幾何誤差の許容値は既存文書に**未定義**であり `accepted=null`。厳密一致は交換契約に整合しているという診断だけで、物理受入ではない。独立生成間の丸め差も今の比較器では拒否し、許容値を測定結果から作らない。法線ノルムの `128 eps` は数値の単位ベクトル表現を検査する丸め規約のみで、幾何差の受入ゲートではない。

- exit 0：`synthetic` の人工データ作成が終了。
- exit 2：`series`、または契約整合の `compare` が終了したが、仕様の受入は未決定。
- exit 1：契約不一致、入力不備、読込/出力エラー。`report.json` の原因を確認する。

本体がこの形式をまだ出していないこと、真の折れ目/段差を現行では拒否することを維持する。交換ファイルに現れない本体内部の隠れた補正まで証明するものではない。独立な元入力、バイナリSHA、元の未加工出力と本体ログをユーザーPCで保持する必要がある。
