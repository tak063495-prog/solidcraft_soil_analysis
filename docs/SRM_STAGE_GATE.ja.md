# 途中工程のFsゲート

未低減の受理状態を起点に、指定した施工工程の終わりで要求率Fs_reqを評価し、後続工程を続けるか止めます。Native SolidCraftの画面には接続していません。一般斜面の正確なFs、連続体安定性、設計安全認定、P6全体の完成を意味しません。

```text
solve-mc-epoch-srm-gate MESH.json GROUND.json CASE.json SPEC.json NEW_REPORT.json
prepare-mc-epoch-srm-gate-job MESH.json GROUND.json CASE.json SPEC.json ACCEPTED_UNITS NEW_JOB
run-mc-epoch-srm-gate-job JOB [ADVANCE_BUDGET]
stop-mc-epoch-srm-gate-job JOB
mc-epoch-srm-gate-job-status JOB
verify-mc-epoch-srm-gate-job JOB
prepare-mc-epoch-srm-gate-recovery PREVIOUS_JOB NEW_JOB
```

```json
{
  "schema": "soil-fem-epoch-srm-stage-gate/1",
  "search": {
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
      "through_stages": 0
    },
    "initial_increment": 0.25,
    "minimum_increment": 0.000244140625,
    "maximum_attempts": 128
  },
  "gate_after_stages": [1, 2]
}
```

`search.search.through_stages` は0です。各ゲートが、その工程数まで完成した未低減Runへ工程数を割り当てます。0以外は計算前に拒否します。`gate_after_stages` は1から全工程数までの厳密増加列です。

各ゲートでは、その工程までの未低減施工を完成させてから[epoch起点の探索](EPOCH_SRM_SEARCH.ja.md)を行います。低減後の材料状態は次の工程へ渡しません。判定は次の三つです。

|効果|条件|後続工程|
|---|---|---|
|`continue`|要求率の結果が `equilibrium_requirement_met` で、停止理由が `factor_tolerance` または `no_physical_upper`|続ける|
|`stop_excluded`|要求率の結果が `independently_excluded`|止める。不合格であり、崩壊認定ではない|
|`stop_indeterminate`|要求率が `indeterminate`、または釣合いを得ていても試行上限・矛盾・未完了の探索|止める。判定不能であり、不合格にも上端にもしない|

試行上限で区間幅が閉じない場合は、要求率で釣合いを得ていても次工程を認可しません。不収束、最小刻み、材料試行上限は判定不能のままです。未低減施工そのものが収束しない場合はゲート判定を作らず、求解エラーとして直前の受理単位を残します。

全ゲートが継続し、残りの工程も終えたときの `decision` は `passed` です。これは要求率の数値的な釣合いを、指定した工程の終わりで得たという意味です。`stopped_excluded` と `stopped_indeterminate` は、止めた工程までの未低減の完成履歴だけを残します。後続工程は開始しません。`nonconvergence_used_as_exclusion_or_collapse` と `continuum_stability_or_design_safety_certified` は常にfalseです。

ジョブの終了コードは、ゲート計算が終わったとき0、受理単位での休止が2、入力・制御・保存の異常が1です。0は合格ではありません。`decision`、各ゲートの `effect`、`requirement`、`stop_reason` を読んでください。保存は施工の受理増分、または探索中の材料試行境界を1単位とします。棄却された材料試行は受理状態を置き換えません。停止はハードリンクの要求を本人が確認し、二つのスロットの完成したポインタだけを復旧します。元ジョブは変更しません。

全ゲートを通過した未低減施工は、同じ入力の連続 `solve-stages` と、計測時間およびその時刻を含むRun全体ハッシュを除いて一致します。途中で止めた接頭は、後続に表示RESETがある場合の材料ID配列を除き、連続計算の対応するフレーム・試行・確定状態と一致します。表示RESETは全Caseの生命周期フラグを変えるため、接頭だけのCaseとは材料IDの保持が異なります。変位・応力・試行の採否は比較します。

対象は既存のepoch探索と同じ乾燥MCと、圧力・境界・掘削・表示RESETの接頭です。水、構造、JOINT、Davis、MCC/IDH、実断面、Nativeの工程画面、一般3D精度、元の2D変位5％基準は未完です。
