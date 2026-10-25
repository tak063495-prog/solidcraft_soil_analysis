# SolidCraft外壁面と平板の対応・限定検証

`plate_surfaces.py` は名前付き外表面セットを、地盤節点に接続した平板へ対応させる。Tri6面は中間節点を含む四つのTri3平板に分ける。元の座標・体積要素・材料・領域・既存の面セットは変更しない。専用の一節点セットで接続するため、異なる土塊の重複座標を座標だけで結び付けない。

面セット、所有ボディ、板断面、初期有効状態、一定/比例面トラクションを `soil-fem-plate-surfaces/1` に指定する。面の外向き・面積を検査し、内部面、重複選択、ID衝突、曲がったTri6、構造節点予算超過を拒否する。既存構造の接続節点が一意なら同じ節点ID・回転を再利用し、同じ地盤節点に異なる構造回転がある場合は自動統合せず拒否する。設置/撤去の `members` に指定した面グループIDを、生成した板ID群へ展開する。

```powershell
.venv/Scripts/python.exe tools/prepare_plate_wall.py runs/my-wall-input --cli tools/solidcraft-0.20.0/solidcraft-0.20.0-windows-x64-portable/solidcraft-cli.exe --mesh-sizes-m 4 3
target/release/soil-fem-3d.exe solve-stages runs/my-wall-input/coarse/ground/mesh.json runs/my-wall-input/coarse/ground/ground.json runs/my-wall-input/coarse/stages.json runs/my-wall-input/coarse/stage-results.json
.venv/Scripts/python.exe tools/export_stages.py runs/my-wall-input/coarse/ground/mesh.json runs/my-wall-input/coarse/stage-results.json runs/my-wall-input/coarse/exports
```

最初のコマンドは、実SolidCraft CLIでボディを作成・STEPへ保存し、粗密メッシュと解決済みの板対応・施工入力を準備する。平衡済み・解析済みとは記録しない。後続の計算には板を含む最新版のソルバーをbuildして使用する。出力先は新しいフォルダー名を指定する。

## 対照モデルと結果

`examples/plate-wall.json`、`ground-plate-wall.json`、`stages-plate-wall.json`、`plate-wall-surfaces.json` を入力にした。地盤は10×10×8 m、γ=20 kN/m³、E=30,000 kPa、ν=0.3。MCの弾性枝に留める強度を明示しており、地盤強度の同定例ではない。底面並進3成分を拘束し、他の面は自由。x=10 mの80 m²の外壁にE=30,000,000 kPa、ν=0.2、厚さ0.2 m、γ=24 kN/m³の板を設置する。自重384 kN、水平面トラクション−5 kPaによる計400 kNを載除荷し、板を撤去する。

記録は `runs/plate-wall-cad-20261026-v2/verification.json`。実行版は `8a3eda114a0d92a46f5854b84e08c010928c513f7b2a0ccf0fc74dfd86da2e7d`、実行ソース555ファイルを固定した。

|項目|粗|細|
|---|---:|---:|
|地盤Tet10 / 節点|156 / 359|335 / 690|
|平板数 / 共有壁面節点|80 / 51|136 / 83|
|板の面積平均水平変位・載荷増分|−1.416447 mm|−1.430086 mm|
|地盤Tri6面で積分した平均水平変位・載荷増分|−1.421187 mm|−1.433552 mm|

粗密の板平均変位差は細側基準で0.95368％。今回固定した2％のペア感度基準は満たしたが、2メッシュ・一つの機能量だけで一般的な空間収束とは判定しない。Drilling係数1e−6/1e−4/0.01を粗側で比較し、全節点変位・全GP応力・全板評価点曲げモーメントの載荷増分L2差は最大9.49e−7で、事前の1％基準を満たした。面法線回りの梁モーメントなど、別の接続条件の感度を保証しない。

4条件の各5フレームを、元入力から別に組み立てた全体3D系と照合した。地盤の全節点・全GP応力/e、構造の全並進/回転・設置基準、全板評価点の膜力・曲げモーメント・エネルギー、反力、力/モーメント、VTU/生値を検査した。粗側では初期化・設置・載荷・除荷・撤去の5地点から別CLIへ再開し、状態・全フレーム・監査・非時間試行履歴が一致した。意図した切戻しを除いた経路とも状態・全フレームが一致した。細側の別CLI再開を今回受入したという意味ではない。

## 継続する要件

共有した全節点の並進は一致するが、Tri6地盤面と四分割したTri3板の補間は節点間では一致しない。平均値の差は粗4.740 μm、細3.466 μmだった。この差を隠して面全体の点ごとの変位連続を受入済みにしない。整合する高次板・弱い界面結合等の定式化と、その収束検証が残る。

これは有限幅の外壁の弾性対照であり、掘削を伴う山留め壁・塑性地盤・水圧/JOINT・内部掘削面の保持側選択・GUI面選択や専用表示の受入ではない。CAD編集後の面IDの永続対応、一般境界/奥行き・メッシュ収束、P1体積投影/3場MCCと部材の接続も継続する。[平板の方式](PLATES.ja.md)と [全残件](P6_STATUS.ja.md)を参照する。
