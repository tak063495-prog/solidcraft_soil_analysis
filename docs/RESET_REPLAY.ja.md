# 表示RESET・状態RESETと独立施工再生

2026-10-24。表示変位のリセット、明示した材料状態への置換、MC強度係数ごとの独立施工再生を実装した。有限変形や、変形後座標へ基準を変更する2D互換RESETとは区別する。

```powershell
cd D:\Work\soil-fem-3d
.\scripts\Invoke-Stages.ps1 -Design examples\construction-fill.json -GroundFile examples\ground-construction-fill.json -StageFile examples\stages-reset-replay.json -ElementOrder 2 -MeshSizeM 4
# 表示された解析フォルダを指定する
.\scripts\Invoke-Replay.ps1 -SourceRun runs\<source-run> -SpecificationFile examples\replay-mc.json
.\scripts\Test-ResetReplay.ps1
```

## 表示だけをリセットする

```json
{"kind":"reset_display","name":"zero_construction_display"}
```

現在の確定変位を`display_reference_displacement_m`に保存する。内部の累積変位、応力、塑性ひずみ、MCCのpc/e、MCの塑性仕事、物理荷重、有効領域は変更しない。材料積分を行わない原子的な操作で、ゼロ増分の帰着によって材料の診断履歴を書き換えることも避けている。

VTUの`displacement_m`と監視値は「内部変位−現在の表示基準」。`total_displacement_from_seed_m`は内部の累積変位、`initialization_displacement_m`は初期平衡化の変位、`display_reference_displacement_m`は現在の表示基準である。初期平衡からのひずみ指標も維持するため、表示RESETで蓄積変形を隠して小ひずみ判定することはない。

RESET後にBIRTHした新規節点は、配置時の位置をその節点の新しい表示基準にする。既存有効地盤と共有する節点の基準は維持する。撤去済み土の表示基準を新しい盛土へ引き継がない。

## 材料状態を明示してリセットする

```json
{"kind":"reset_state","name":"new_initial_state","bodies":["Ground"],
 "states":{"Ground":{
   "effective_stress_compression_positive_kpa":[0,0,0,0,0,0],
   "plastic_strain_compression_positive":[0,0,0,0,0,0],
   "void_ratio":0.7,"equivalent_plastic_strain":0,
   "plastic_work_kpa":0,"last_return":"elastic"
 }}}
```

上の状態はMCの例。選んだボディごとの全材料状態を入力し、現在材料の型・降伏・pc/e整合を検査する。MCCは正の平均有効応力・pcと整合する間隙比を指定する必要があり、上のゼロ応力MC例をMCCへ自動転用しない。

材料番号、基準座標、内部変位、表示基準を維持して全積分点の材料状態を置き換える。新しい種状態の内力と実荷重の不平衡を仮保持荷重で受け、同じ非線形FEMで段階的に釣合いを取り直す。自重下で応力をゼロに置き換えると、その荷重で再び変形する。表示だけをゼロにしたいときにこの操作を使わない。

MATSETと同様に、置換直後・平衡増分前の状態、仮保持荷重、λ=0の残差を`transition_audits`に保存する。切戻しでは置換前の確定状態へ戻り、再試行で改めて同じ状態を作る。表示RESETはそれ自身の確定境界、状態RESETは受理荷重増分の途中から保存再開できる。

両RESETを使う結果は`soil-fem-stage-result/6`。既存のRESETなし入力では従来の結果スキーマ・数値を維持する。

## MC強度係数ごとの独立再生

```json
{"schema":"soil-fem-replay/1","material_ids":[1,2,3],"factors":[1,2,3],
 "psi_policy":"cap","tension_policy":"keep"}
```

各係数Fで元入力を複製し、選択材料の`c_F=c/F`、`tanφ_F=tanφ/F`を設定する。**初期平衡前の元入力から、その係数で初期平衡と全施工を再計算する**。前の係数の変位・塑性履歴・保持荷重・切戻し状態は持ち越さない。MATSET用に宣言した将来材料も、対象番号へ含めるか明示する。

ψは`keep`、`zero`、`cap=min(ψ,φ_F)`から選ぶ。関連流れのψ=φ条件などに反する指定は入力棄却として記録する。引張カットオフは`keep`または`divide_by_factor`を選ぶ。Davisでは入力c/φとψ方針を適用した後、Davisの有効強度を再評価する。E/ν・重量・初期水条件を強度係数で変更しない。

現在は全材料がMCの地盤を対象にし、MCCを含む指定は拒否する。係数は1〜100、重複なし、最大50試行。係数1で元のψ/引張方針を保つ場合は通常解析と完全一致する。`zero`方針が元の非ゼロψを変更する場合は、係数1でも別の流れ則になる。

`Invoke-Replay.ps1`は受入済みの元解析パッケージのハッシュを検査し、入力・元結果・仕様と実行exeをジョブへ固定する。元結果の形成履歴を新しい係数へコピーする処理ではない。出力は`replay-results.json`と`trials/trial-NNN/`。完了試行だけに独立した力・モーメント・履歴検査を行い、受入済みVTU/PVDを作る。

生のCLIは`replay-stages MESH.json GROUND.json CASE.json SPEC.json NEW_REPORT.json`。未完了・材料エラー・不正方針はそれぞれ区別して保存し、CLIは終了コード1を返す。数値未完了で得た最後の確定状態は診断用に保持する。不正なパラメータ方針では実行されていない調整材料や解析結果を作らない。

## 施工工程の終点を指定する

```text
soil-fem-3d replay-through-stage MESH.json GROUND.json CASE.json SPEC.json THROUGH_STAGES NEW_REPORT.json
```

THROUGH_STAGESは元CASEの施工工程数までの整数で、0は初期平衡のみ、1は最初の施工、最大値は全施工です。元の完全なCASEとそのidentity、要求工程数、実際に再生した部分CASEをsoil-fem-replay-stage-gate/1へ区別して保存します。除外した後半工程の診断用棄却注入は実行せず、元CASEには保持します。

各係数は低減済みの材料で初期平衡前から独立再生します。無低減の施工状態を一度作り、そこから強度だけを落とす機能ではありません。仕様は既存soil-fem-replay/1で、MCC等は拒否します。不正な工程数は計算前に拒否し、既存出力を上書きしません。全試行が終点へ到達した場合だけ終了0とし、数値未完了・入力棄却・材料/求解エラーは区別して報告します。

終点到達は安全率の認定ではありません。stability_classificationはnot_evaluatedのままで、Fs探索・Fs_reqの途中施工ゲートは後続要件です。選択した工程終点と、SRMの物理的合否判定を混同しないでください。

## SRMの残り

この再生機能はSRMの基盤で、Fs探索・Fs_reqの施工途中ゲート・物理的上端候補の判定はまだ行わない。出力の`stability_classification`は`not_evaluated`。完了はFEMの収束を意味し、安全の認定ではない。非収束・線形求解失敗・入力棄却から崩壊やFs上端を生成しない。

確定した途中起点からの再生選択、2D版の変形後基準座標RESETとの互換、Fs探索と独立な破壊基準・粗密/奥行き検証は後続要件として保持する。[検証記録](RESET_REPLAY_VERIFICATION.ja.md)／[P6一覧](P6_STATUS.ja.md)
