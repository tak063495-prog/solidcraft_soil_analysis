# SolidCraftの局所掘削と3D地盤変位

2026-10-23。水平多層地盤の初期状態を保ったまま、SolidCraftで指定した閉じた掘削用ソリッドとの交差領域だけを撤去する。RustのMCC全体非線形FEMで、3方向変位と応力・pc・間隙比・塑性ひずみを計算する。

## モデル作成と実行

1. SolidCraftで掘削前の地層ボディを作る。標準例は10×10×8m、上下各4mの2層。
2. 地層を直接切り取らず、独立した閉ソリッド`Pit`を重ねて配置する。標準例はX/Y=3〜7m、Z=6〜9mの箱。地盤との交差は4×4×2m=32m³。地表より上に出た1mは解析領域に加えない。
3. `ground-local.json`の`layers`に材料地層名、`excavation_bodies`に掘削用ボディ名を登録する。名前は重複不可。未登録ボディがあるCADは拒否する。
4. `stages-local-sparse.json`の`kind: "excavate"`、`bodies: ["Pit"]`で工程を指定する。局所掘削スクリプトはこのCSR用入力を既定とする。`stages-local.json`は密LU比較用として保持する。

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-LocalExcavation.ps1
.\scripts\Invoke-LocalExcavation.ps1 -Design .\my-model.json -GroundFile .\my-ground.json -StageFile .\my-stages.json
.\scripts\Invoke-LocalExcavation.ps1 -Live -Connect 127.0.0.1:7879
.\scripts\Test-LocalExcavation.ps1
.\scripts\Test-LocalExcavationAcceptance.ps1
```

`-Live`では起動中SolidCraftのモデルをスナップショット化する。CADの編集内容・材料値・工程は実行時に読み取る。解析ボタンをSolidCraftのUIへ組み込む改変はしていない。

現行Tet10の箱・円柱掘削、独立3D比較、塑性履歴を伴う切戻しと別プロセス再開を一括で確認するコマンドは`Test-LocalExcavationAcceptance.ps1`。入力・exe・検証ログを実行別フォルダに固定し、合格時に`final-validation.json`とソースのZIPを保存する。[現行受入記録](LOCAL_EXCAVATION_ACCEPTANCE.ja.md)を参照。

## CADから撤去要素への対応

`cad_layers.py`は各地層と各ツールをボディ名でSTEP出力し、SHA256・CAD体積・境界箱を`bodies.json`へ保存する。出力前後の設計一致も確認する。

`mesh_layers.py`は全地層とツールをOpenCASCADEで同時fragmentし、Boolean演算の祖先情報から各子体積の元地層とツールへの所属を取得する。境界箱・要素重心だけで掘削対象を近似しない。地層同士の正体積重複は拒否し、ツールと地層の重複は撤去要素集合へ変換する。

元地層の`region_ids/material_ids`は変更せず、`excavation_sets`に要素番号を格納する。ツールのみの子体積は捨てる。CADと各地層のBoolean後体積が一致し、切断面の節点を共有する共形メッシュを作る。局所掘削の既定は曲面中間節点を保持するTet10。`-ElementOrder 1`で従来のTet4を使用できる。どちらも曲面積・撤去体積に離散化誤差がある。ツールとの交差CAD体積は`mesh.provenance.json`にも記録する。Tet10の積分則・品質検査・体積計測方法は[Tet10仕様](TET10_BOUNDARY.ja.md)を参照。

撤去領域の各連結成分は、現在露出している地表または既存掘削面につながっていなければならない。残存地盤は面を介して支持底面につながっていることを検査する。埋没空洞、既に撤去された領域との重複、全地盤撤去は拒否する。複数工程は重複しない別ツールで指定できる。

## 応力解放と荷重

確定済み応力から撤去要素の内力を組み立て、仮の保持荷重を`H=F_old−F_new−F_int,removed`とする。各掘削増分で`F(λ)=F_new+(1−λ)H`を解く。λ=0は従前の残存地盤とつり合い、λ=1で仮保持力を完全に消す。撤去土の自重はF_newに含めない。

`excavation_audits`に撤去要素数・体積・初期保持荷重・λ=0の自由度残差を保存する。解放荷重は撤去面の応力を等価節点力として解放するもので、切断直後から残存地盤の応力を0へ書き換える処理ではない。乾燥した自由掘削面は弱形式の自然境界条件として外力0となる。有限要素の点ごとの面応力が厳密に0になることを保証する条件ではない。

局所掘削には`pressure_on_original_surface_only: true`を必須とする。`remaining_surface_pressure_kpa`は元の地表Z高さに残った面だけへ作用し、新しい床面には上載圧を載せない。従来の全層ベンチマークはこのフラグが省略されるため、露出した水平面に指定圧を作用させる以前の仕様を保つ。

標準例は乾燥地盤。水位を指定して初期水域の内部を掘削すると、その水位まで湛水した掘削を表す。新たな床面・側壁に外水圧を作用させ、土内部には間隙水圧を加える。現行版は曲面Tri6の部分水没水圧積分と、指定水位の載除荷工程にも対応する。[静水圧仕様](HYDROSTATIC.ja.md)を参照する。K0初期化の水位・地層境界の整合条件は維持する。未知の自由水面、掘削排水の浸透・過剰間隙水圧・圧密の時間変化は本受入の範囲外。

## 切戻しと結果

Newtonの試行状態は確定履歴から独立して作り、失敗した試行の変位・応力・pc・間隙比・塑性ひずみ・有効要素・荷重・解放進行度を確定しない。刻みを縮小して同じ確定状態からやり直す。撤去済み要素の材料履歴は凍結する。刻み下限で失敗した場合は`complete=false`と`last_committed`を保存し、合格結果のVTUを出力しない。

`Pit:floor`は残存地盤から掘削領域へ向く三角形のうち法線のZ成分がノルムの半分を超える面、`Pit:wall`はその他の内部切断面、`Pit:surrounding_surface`は未撤去の元地表面から作る節点集合。傾斜面の監視位置を厳密に指定したい場合は`node_sets`を別途与える。複数ツールの集合はそれぞれの単独撤去形状から作られ、集計時には現工程の有効節点で絞る。

`monitor_sets`に監視名と節点集合名を指定する。平均変位は節点の算術平均で、面積平均ではない。最小/最大変位もJSON/CSVに保存する。Z上向き正、沈下は負、隆起は正。周辺地表に沈下・隆起が混在する場合は平均値だけで判断しない。工程増分とK0からの累積変位も区別する。

`results/stages.pvd`をParaViewで開くと、撤去領域を除いた3D形状上で`displacement_m`のZ成分、Warp By Vector、材料点の応力・硬化状態を表示できる。`stage-001.vtu`が標準1工程の掘削後結果。`stages.csv`は監視結果、`attempts.csv`は切戻し記録。工程番号は擬似時間であり実時間ではない。

独立比較の際には、積分を細かくしたRust計算も`reference-refined-result.json`と`reference-refined-exports/`へ別保存する。比較報告のRust数値はこちらに対応する。通常実行の`stage-results.json`、`results/`は比較のために上書きしない。

## 検証の範囲

`Test-LocalExcavation.ps1`はRust全試験、Gmsh実体積分割5試験、SolidCraftからの箱形局所掘削のTet4回帰、2m/1.5mメッシュでの独立3D比較、塑性MCCの独立1D回帰計算を実行する。曲面Tet10・全積分点履歴・一般境界条件・別プロセス再開は`Test-Tet10AndBoundary.ps1`で検証する。

独立3D比較は、深さからK0初期応力を別に作り、解析的なMCC弾性枝の速度剛性と陽的中点積分で全体FEMを解く。Rustの数値接線・局所return mapping・Newton変位を参照解の生成に使わない。同じ空間離散化の全有効節点3方向変位、全有効要素の応力と間隙比を比較し、参照積分100/200分割とRustの積分細分化を確認する。参照計算が降伏面に接近した場合は拒否するため、局所掘削の塑性域まで独立3D計算が検証済みとはしない。

塑性を伴う局所掘削の履歴復元は、プレロード→除荷→局所掘削で試行収束後に棄却を挿入し、同じ確定刻み列の正常計算と最終確定履歴が完全一致することで検証する。構造化288要素の3D地盤を使い、材料履歴・有効要素集合・載荷状態が実際に変わった試行を棄却する。実際のNewton失敗で刻み下限に達した際の初期状態保持、水中掘削の床面・側壁水圧も試験する。

水平直方体層の初期地盤、静的・小ひずみ、等方MCC、固定水圧という適用範囲は維持する。CSR疎求解と保存再開を追加し、11480要素の局所掘削まで検証を拡張した。[疎求解・再開仕様](SPARSE_RESTART.ja.md)を参照。閉ソリッドによる切断形状への対応は一般斜面のK0初期化対応を意味しない。支保や接触、せん断帯正則化、掘削底面の破壊・安全率、現場定数の同定は別途必要。[MCC定式化の参考資料](https://docs.itascacg.com/itasca920/common/models/camclay/doc/modelcamclay.html)。

空間離散化の感度は`spatial-sensitivity.json`へ別出力する。同一メッシュ上の独立計算との一致から、メッシュ収束や実地盤に対する精度を推定しない。[今回の実行記録](LOCAL_EXCAVATION_VERIFICATION.ja.md)に、2段階メッシュの差と適用限界を示す。
