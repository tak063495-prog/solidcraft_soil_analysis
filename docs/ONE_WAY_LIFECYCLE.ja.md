# 水理領域の再配置・材料変更と履歴

`one_way_construction.rs`へ、初期非有効ボディ、物理時刻0の掘削、湿潤MCの無応力BIRTH、MATSETを追加した。既存の局所DEATHと全体Newton・Richardsを用い、施工解放のλと物理時刻を分ける。各試行は材料・水理・領域・水量台帳をまとめて確定する。小ひずみ、固定水理間隙率と乾燥かさ重量を使うONE_WAYである。

## 操作

局所掘削はSolidCraftの別ソリッド`Pit`を使い、Booleanの実際の交差要素を撤去する。[局所掘削の入力](LOCAL_EXCAVATION_CURRENT.ja.md)と[水理DEATH](ONE_WAY_CONSTRUCTION.ja.md)を参照。新しい実行ラッパーはSolidCraftのDesign/Liveまたは準備済みメッシュを受け取る。

```powershell
pwsh -NoProfile -File scripts/Invoke-OneWayConstruction.ps1 `
  -Design examples/local-excavation-cylinder.json `
  -GroundFile examples/ground-local.json `
  -CaseFile YOUR_MATCHING_3D_CASE.json -ElementOrder 2 -MeshSizeM 4
```

Caseの期間・水理/外水境界・材料番号は地形に合わせて明示する。柱用`floor`を任意CADの面名とみなさない。メッシャーは二つの地層の共有面を`interface:LowerClay:UpperClay`のような辞書順のボディ名で保存する。初期非有効の上層がある場合、この共有面を初期の地表境界として指定できる。現在のボディ名での対応で、CAD編集後の永続参照を保証しない。

最小柱例は次で実行できる。

```powershell
target/release/soil-fem-3d.exe solve-one-way-construction `
  examples/seepage-column-tet10-mesh.json examples/ground-one-way-lifecycle.json `
  examples/one-way-lifecycle.json runs/NEW-lifecycle.json
pwsh -NoProfile -File scripts/Test-OneWayLifecycle.ps1
```

ラッパーの`-PreparedModel DIR`はDIRの`mesh.json`を使う。実行ファイルを出力先へ固定コピーし、解析後にVTU/PVDを再読込して確認する。元入力、固定実行版、最後の確定状態を保存する。出力先は毎回新規とする。

## 初期非有効、BIRTHと材料変更

外側入力へ`initially_inactive_bodies:["UpperClay"]`を指定すると、初期平衡と水理ICを有効領域だけに構成する。後の期間の水頭を初期状態へ流用しない。初期水頭から作る力学平衡と、正時間の水理平衡は別である。

`transitions`へ次の形式のイベントを置く。各期間の開始で掘削または一つのtransitionを実行できる。時刻0でも偽の微小dtは挿入しない。

```json
{
  "name":"place", "before_period":"born",
  "water_boundary_policy":"explicit_period_boundaries",
  "surface_traction_policy":"explicit_period_reservoirs",
  "action":{
    "kind":"birth", "bodies":["UpperClay"],
    "states":{"UpperClay":{
      "effective_stress_compression_positive_kpa":[0,0,0,0,0,0],
      "plastic_strain_compression_positive":[0,0,0,0,0,0],
      "void_ratio":0.7, "equivalent_plastic_strain":0,
      "plastic_work_kpa":0, "last_return":"elastic"
    }},
    "placement":"cad_new_nodes", "initial_head":{"head_m":0}
  }
}
```

BIRTHは非有効ボディと新しいゼロ応力MC状態を必要とする。`cad_new_nodes`は新たに使う節点の変位を0にし、既存共有節点の変位を保つ。`current_displacement`は保存された節点変位を配置基準にする。配置準備は最初の試行だけに適用し、解放の各増分で繰り返さない。MCCはp=0で体積剛性が0になるため、無応力BIRTHを受け付けない。

MATSETのactionは`kind:"material_set"`, `bodies`, `material_id`, `history`, `water`を持つ。材料番号を地盤の`construction_material_ids`に登録し、対応する水理材料則も指定する。`history:{"mode":"preserve"}`は応力と履歴を保持する。`mode:"replace",states:{...}`はボディごとの明示状態へ置き換え、適合性を検査する。

水理方針は次のいずれかを明示する。

|water.kind|作用|
|---|---|
|`preserve_head`|旧Hを保持し、新しい保水則で水量を求める。差を搬入/搬出水量として記録|
|`preserve_water`|セルの実水量を保持し、新しい保水則を積分してHを逆算。`reference_head`を必須とする|
|`initialize_head`|`value`の明示Hで新しい水理状態と搬入/搬出量を作る|
|`inactive`|全選択セルが非有効のMATSET専用。水頭や水量を作らない|

Ss=0の飽和域はHが水量から一意に決まらない。`preserve_water`は与えた参照Hがその水量に適合する場合のみ採用する。残留含水量以下など達成できない水量を要求すると、その工程起点で失敗し、直前の確定状態とアーカイブを保持する。将来の実状態を仮造して適合性を先判定しない。

## 水量、応力解放、保存

実水量は`theta_s V + W`で、Ssによる正圧貯留も含む。累積監査は

`残存水量 + 搬出水量 − 搬入水量 − 初期水量 − 積分した正味流入`

で行う。変更していないセルのH/W/局所Q、既存面のトレースは完全に移す。新規面だけ明示ICからトレースを作る。変更セルのQは新しい材料/面作用素により再構成する。この瞬間を水理平衡とは呼ばず、その後の正dtで解く。

準備後の状態と新しい外力から仮保持力を組み立て、λ=0の平衡を監査してから0へ解放する。履歴を保持したMATSETは旧応力を引き継ぎ、新しい剛性による増分を解く。BIRTHが受け持つ水重量、χuの内力、外水圧を別々に組み立てる。棄却時は水理材料ID、間隙率参照、配置基準、塑性履歴、仮保持力、台帳を全て戻す。

撤去GPは撤去時の材料IDと全状態を別保存する。非有効域の材料変更や再配置が後続しても、旧MCC応力・pc・e・塑性履歴の記録を上書きしない。材料変更前の履歴も別VTUへ出す。保存再開は入力/版照合と経路再生を伴い、解放途中の保存から同じ次増分へ続ける。

## 検証範囲

Rustの追加12試験はTet4/10の飽和柱BIRTH理論、材料水量移送、時刻0のMCC掘削、正圧貯留、実際に塑性履歴が変わった試行の棄却、保存改変拒否、非有効MATSET→BIRTHによる旧MCC記録の保持を含む。総合スクリプトは全Rust試験と従来の施工/固定域/水理単独/既存結果検査を同じ凍結ソースと実行版で取り直す。

独立側は保水積分・移送水量・未変更主状態・トレース・旧材料履歴、荷重と材料内力、局所RT0収支を再組立する。別の増分弾性FEMは入力のK0と自分の変位・応力を保持し、BIRTH配置/材料変更後の剛性と履歴を更新する。独立積分5/7/9を比較し、変位増分のノルムで誤差を評価する。水頭はRustの入力受渡し値を使うため、独立Richards PDEの検証ではない。一般3D塑性変位の独立解でもない。

湿潤JOINT/部材、無応力MCC誕生、連成圧密、湛水動態、一般3D塑性/有限幅過渡前線の独立比較、SolidCraft内UIなどP6全体の残作業は継続する。

## 限定受入結果

2026-10-24の[総合受入](../runs/lifecycle-acceptance-20261024-v2/final-validation.json)と[integrity.json](../runs/lifecycle-acceptance-20261024-v2/integrity.json)。実行版SHA256は`46a2ec1895e8f3991047f5f7abad9dba977f2335d7206021fa1b563a2e55c6ce`。開始前の385ソースと実行版を固定し、2,589証跡ファイルを束ねた。

- Rust **197試験**、CAD分割5件、接合面分割2件、整形・Clippy警告拒否・releaseビルドが合格。
- 新施工**33地点**の別プロセス再開で全状態・全フレーム・試行・アーカイブが完全一致。同じ版で従来施工28地点、固定域連携29地点、水理単独8地点と既存6ケース完全一致も取り直した。
- [実CADの追加2地点](../runs/lifecycle-cad-checkpoint-20261024-v1/verification.json)はBIRTH/MATSET途中の非ゼロ保持力と固定時刻を保存し、完全一致で再開した。[材料変更後の掘削1地点](../runs/lifecycle-cut-after-matset-20261024-v1/verification.json)でも、現在の材料3の保水則・撤去水量・応力履歴を独立再組立し、途中再開の完全一致を確認した。この別票は`initialize_head`方針も含む。合計36地点は33地点と追加3地点を分けて集計する。
- 塑性履歴が変化した候補の棄却と、正常な受理刻み列との一致を確認。不可能な水量移送は実際の工程起点で停止し、直前の確定状態と旧MCCアーカイブを保持した。
- 飽和柱のBIRTH沈下は理論値4.691176mmに対し、Tet4/Tet10とも絶対差5×10⁻¹⁴m未満。正圧の比貯留も搬入水量へ反映した。
- ラッパーからSolidCraft CLIを再実行し、595節点・312 Tet10を生成。位相・集合・材料対応は完全一致し、再取得時の座標差は最大7.14×10⁻¹³m（許容10⁻¹¹m）だった。結果・生GP・旧履歴のVTU/PVD再読込も合格。

独立増分弾性FEMによる、実CADの全有効節点3成分のイベント増分比較。各増分自身のノルムで正規化する。ほぼ変位しない物理時間ステップをこの表に混ぜない。

|イベント|増分の相対L2差|独立積分7→9の増分差（比率）|
|---|---:|---:|
|湿潤MC BIRTH|0.00011305％|8.96×10⁻¹¹|
|履歴保持MATSET|0.00017633％|5.19×10⁻¹⁰|

水頭は受渡し値であり、独立水理PDEや一般MCC塑性3Dの全面照合ではない。12ケースの独立再組立の最大力残差は6.93×10⁻⁹、水量台帳誤差は最大3.13×10⁻¹³m³だった。この集計は独立再組立側で、解析が保存した収支と区別する。

[再配置・材料変更の3D変位増分と水量の図](../runs/lifecycle-figure-20261024-46a2ec18-v2/lifecycle.png)は原形状を着色し、全有効節点の生の最小・最大値を色範囲に含める。新規節点を0にする配置基準を両計算へ適用した。作図ソース・入力ハッシュ・目視確認も保存する。局所Pit撤去の図とは別の施工例である。

検証後の文書と図、追加3地点は[現行ソース照合票](../runs/lifecycle-current-source-binding-46a2ec18.json)に対応付ける。P6全体の完了とは判定していない。

上記のソース照合は46a2ec18固定版に対するもので、後続の[湿潤JOINT](JOINT_WATER.ja.md)追加ソースを含めない。
