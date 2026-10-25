# 盛土の無応力誕生と材料変更

2026-10-24。一般初期平衡→BIRTH→MATSET→載除荷→DEATH→材料交換→再盛土を全体FEMへ接続した。`birth`と`material_set`は、MCC/MCの既存非線形Newton・切戻し・チェックポイントを使う。小ひずみ・固定メッシュの施工解析である。

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-Stages.ps1 -Design examples\construction-fill.json -GroundFile examples\ground-construction-fill.json -StageFile examples\stages-construction-fill.json -ElementOrder 2 -MeshSizeM 4
.\scripts\Test-Construction.ps1
```

SolidCraftでは、初期地盤と、後で盛る土を別の閉じた材料ボディとして作る。両方を共形メッシュに含め、`initially_inactive_bodies`で将来の盛土を初期状態から除外する。撤去用ツールと盛土材料ボディは役割が異なる。盛土の外形を掘削ツールとして指定しない。

## 初期状態

`initially_inactive_bodies: ["UpperClay"]`を施工ケースに指定する。最初から無効な領域を持つ解析は`initialization.method: equilibrate`を要求する。将来領域の形状・材料・種履歴も入力し、初期平衡の自重、内力、自由度からは除外する。全メッシュのCAD体積・点数・材料対応の検査は維持する。

既存地盤は底面支持へ面でつながる必要がある。現在の連結検査は水平モデル底面を使う。初期化中の全GP履歴、形成した初期応力、施工変位の基準は従来どおり保存する。

## BIRTH

```json
{"kind":"birth","name":"place_fill","placement":"cad_new_nodes","bodies":["UpperClay"],"states":{"UpperClay":{
  "effective_stress_compression_positive_kpa":[0,0,0,0,0,0],
  "plastic_strain_compression_positive":[0,0,0,0,0,0],"void_ratio":0.7,
  "equivalent_plastic_strain":0,"plastic_work_kpa":0,"last_return":"elastic"
}}}
```

選んだボディの全要素が無効であることを要求する。`placement`を明示し、`cad_new_nodes`では既存有効地盤の共有節点の変位を保持し、新たに有効になる節点をCAD基準位置へ置く。これにより再盛土へ撤去済み土の古い節点変位を持ち込まない。`current_displacement`では無効節点の保管変位も含む現在位置へ置く。どちらもその配置時点を無応力の起点にし、過去の変位を新しい材料のひずみとして加えない。

次に有効自由度、接線、自重、表面荷重を再構成し、盛土の重量を段階的に作用させる。選択体が底面につながる有効地盤から離れている場合、部分的に有効な体へ重複誕生を指定した場合、古い塑性履歴を無応力誕生へ混ぜた場合は拒否する。掘削後の再盛土も同じ規則で新しい状態を明示する。

現在の無応力BIRTHは乾燥したMC材料に対応する。MCCは平均有効応力p=0で体積剛性Kもゼロになり、現定式化では無応力状態から開始できないため拒否する。圧力を小さな値へ切り上げて代用しない。水中の無応力誕生も、初期有効応力・外水圧・排水条件を別途定義するまで拒否する。既存地盤にMCC、盛土にMCを使う混在解析は可能。

## MATSET

```json
{"kind":"material_set","name":"improve_fill","bodies":["UpperClay"],
 "material_id":3,"history":{"mode":"preserve"}}
```

将来使う材料も`ground.materials`へ定義し、初期ボディに対応しない番号は`construction_material_ids: [3]`で宣言する。施工中の各要素の材料番号を状態として保存し、構成則・単位体積重量・降伏検査・出力へ反映する。メッシュの元の材料番号は初期CAD対応として保持する。

`preserve`はそのままの応力・pc・間隙比・塑性履歴を新しい構成則へ渡す。新材料の型・降伏条件・MCC基準線との整合を検査し、不整合なら拒否する。モデルの型をまたいだ自動履歴変換は行わない。Eを変更しただけでは過去の変位から応力を作り直さず、以後の増分から新しい剛性を使う。

`history: {"mode":"replace","states":{...}}`では、選んだ各ボディの新しい全材料状態を明示する。有効要素と無効要素の両方に指定でき、無効な撤去体を埋戻し材料へ交換する用途も扱う。新しい状態も構成則の整合検査を受ける。これは明示的な履歴置換であり、履歴継承とは区別して監査記録に残す。

密度を変更したときは現在材料の自重を再構成する。応力や密度の変更による不平衡は、仮保持荷重でλ=0の平衡を作り、λ=1で物理荷重だけになるまで解放する。試行中に変更する材料番号、履歴、有効要素、荷重はまとめて扱い、不受理時には確定状態へ戻す。

## 保存・出力

施工機能を使う結果は`soil-fem-stage-result/5`。確定状態の`element_material_ids`に現在の材料番号、`transition_release_load_kn`にBIRTH/MATSETの解放途中荷重を保存する。初期化の仮保持荷重、掘削の応力解放荷重とは区別する。

`transition_audits`には選択要素、材料変更/誕生直後・平衡増分前の全状態、仮保持荷重、λ=0の残差を保存する。チェックポイントでは元工程状態からこの監査を再構成して照合する。BIRTH途中、MATSET途中、撤去・再盛土途中から連続計算と同じ履歴で再開できる。

VTUのセル`material_id`はその工程の現在材料。全積分点の生の材料履歴と体積平均値は別ファイル/フィールドで出力する。積分点VTUにも現在の材料番号を保存する。内部の累積変位を保持し、監視・表示は初期平衡からの施工変位を使う。BIRTHの工程変位増分は配置後の無応力起点から測り、撤去済み節点の古い変位から差し引かない。その起点は監査の`prepared_state.displacement_m`で確認できる。固定メッシュなので、盛土の体積と重量は基準形状の値を使う。

## 残る施工機能

表示だけのRESET、明示状態置換/再平衡、原入力からのMC施工再生を追加した。[RESET・再生仕様](RESET_REPLAY.ja.md)を参照。2D版のRESET_U/RESET_STRESSは応力・履歴も消し、後者は基準座標も変えるため、3Dの表示リセットと同一視しない。2D互換・基準座標更新、途中起点の選択、SRMのFs探索/ゲート、構造・接触・水理の誕生/除去は未実装。P3a全体、P6全体の完了とは区別する。

[検証記録](CONSTRUCTION_VERIFICATION.ja.md)／[初期平衡](GEOSTATIC.ja.md)／[P6の残要件](P6_STATUS.ja.md)
