# 未低減の施工終了履歴から始める乾燥MCのSRM探索

本体へ接続した経路です。既存の新規初期化からの探索 `search-mc-srm` と、その保存ジョブは変更していません。一般斜面の正確なFs、連続体安定性、設計安全認定、Native UI、P6全体の完成を意味しません。

```text
search-mc-epoch-srm MESH.json ORIGINAL_RUN.json SPEC.json NEW_REPORT.json
prepare-mc-epoch-srm-job MESH.json ORIGINAL_RUN.json SPEC.json MATERIAL_ATTEMPTS_PER_SAVE NEW_JOB
run-mc-epoch-srm-job JOB [ADVANCE_BUDGET]
stop-mc-epoch-srm-job JOB
mc-epoch-srm-job-status JOB
verify-mc-epoch-srm-job JOB
prepare-mc-epoch-srm-recovery PREVIOUS_JOB NEW_JOB
```

`ORIGINAL_RUN.json` は、同じ実行ファイルの `solve-stages` で完成した未低減の施工Runです。起点のメッシュ・入力・実行ファイルID・確定状態と、再構成した物理外力・活動要素・ゼロ変位拘束を検査します。異なる実行ファイルで作られたRunを、自動で現在版へ読み替えません。

各率は同じ元Runから新しい材料経路を作り、c/F、atan(tanφ/F)、明示したψ/引張上限の方針へ材料を段階的に変更します。直前の率の塑性履歴を次の率へ使いません。荷重・施工段階・活動領域・表示RESETの基準は元Runから保持し、撤去済み積分点の履歴も凍結します。元のGroundを新材料へ書き換えた通常Runは作りません。

```json
{
  "schema": "soil-fem-epoch-srm-search/1",
  "search": {
    "schema": "soil-fem-srm-search/1",
    "material_ids": [1],
    "psi_policy": "zero",
    "tension_policy": "keep",
    "maximum_factor": 2.0,
    "required_factor": 1.1,
    "absolute_factor_tolerance": 0.002,
    "maximum_trials": 20,
    "through_stages": 1
  },
  "initial_increment": 0.25,
  "minimum_increment": 0.000244140625,
  "maximum_attempts": 128
}
```

要求率は利用目的に合わせて指定し、規定値や設計安全率として自動採用しません。`through_stages` は与えた完成Runの全工程数と一致させます。途中工程を調べる場合は、その工程までを同じ実行版で完成させたRunを用意するか、[工程ゲート](SRM_STAGE_GATE.ja.md)を使います。

分布拘束では `/2` の `virtual_displacements_m` を指定できます。既存と同じ全積分点の仮想仕事・平均応力の必要条件を使い、独立に排除できた率だけを上端へ使います。候補が不適用、Newtonや線形反復が不収束、材料経路の最小刻みや試行上限へ到達した場合は判定不能です。不収束を破壊側の上端へ変換しません。

出力 `soil-fem-epoch-srm-search-result/1` は、元の未低減Runを一度保持し、各試行を `trials[i].evidence` と `trials[i].epoch` に分けます。`evidence.run` は常にnullです。完成した材料経路だけを現材料の力・降伏・変位拘束で再検査して釣合い成立とします。未完の経路に残る最後の受理状態を、完成した率の釣合いへ読み替えません。

`requirement=equilibrium_requirement_met` は要求率の数値的な釣合いを得たという限定した意味です。`independently_excluded` と `indeterminate` を区別します。区間の下側は釣合いを得た率、上側は独立に排除した率です。矛盾する釣合い/排除の証拠がある場合は判定不能とし、区間を返しません。

保存ジョブは材料試行の境界で二つのスロットへ交互に書き、同期後にポインタを公開します。停止要求はハードリンクで公開し、現在の材料試行を終えて受理済み状態を残します。棄却された試行の応力・履歴を受理状態へ移しません。終了コードは、探索計算の完了が0、受理境界での休止が2、入力・制御・保存の異常が1です。計算の完了は安全率の達成ではありません。復旧は元ジョブを変更せず、新しいディレクトリーを作ります。

現対象は既存と同じ乾燥MC、明示的なゼロ変位拘束、対応済み施工接頭部です。水・構造部材・JOINT・Davisの単調探索・MCC/IDHからの変換・材料置換などは拒否します。材料経路のλは補間座標であり、安全率ではありません。低減後の専用状態を通常の施工保存点へ変換しません。

検証は単軸の理論比4/3、要求率の実評価、未完を未知のまま保持、塑性プレロード・除荷後の柱、工程起点と各率の履歴分離、撤去済み点の凍結、協調停止・プロセス終了後の復旧と連続探索との報告一致、不正起点と封印改変の拒否です。比較で計測時間以外の力学履歴を見ます。SolidCraftの画面、実断面、一般3D精度、元の2D変位5％基準、全P6は別途未完です。
