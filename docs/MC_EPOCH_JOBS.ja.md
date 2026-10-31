# 施工後のMC材料継続ジョブ

完成した未低減の施工結果を起点に、材料定数だけを段階的に変更します。荷重・境界変位・活動要素・施工管理情報を固定し、各試行で全体釣合いを解きます。受理した材料定数と履歴を同時に確定し、材料・求解・線形計算の失敗では両方を戻して刻みを半減します。

対象は乾燥した全MC地盤です。構造部材・接合・水理を含む継続、MCC/IDHの材料変換は未対応です。通常の地盤初期化・施工保存点の検査は保持します。

## CLI

```text
prepare-mc-epoch-job MESH.json ORIGINAL_RUN.json SPEC.json NEW_JOB
run-mc-epoch-job JOB [ATTEMPT_BUDGET]
stop-mc-epoch-job JOB
mc-epoch-job-status JOB
verify-mc-epoch-job JOB
prepare-mc-epoch-recovery PREVIOUS_JOB NEW_JOB
```

最初に同じ実行ファイルの`solve-stages`で未低減の施工Runを完成させます。`prepare-mc-epoch-job`は入力と実行ファイルを新しいジョブへ固定し、最初の再開点を保存します。以後は`NEW_JOB/solver/soil-fem-3d.exe`を使います。

`run-mc-epoch-job`の終了コードは、完了が0、中断保存が2、数値未完・入力/実行エラーが1です。`report.json`の`complete`と`stop_reason`も確認してください。`stability_classification`は常に`not_evaluated`です。

`ATTEMPT_BUDGET`はこの実行で行う最大試行数で、受理・棄却の両方を数えます。指定回数後に中断保存するため、段階的な確認にも使えます。完了や数値未完が先に確定した場合は、その終了状態を返します。

中断要求は別の端末から`stop-mc-epoch-job`を実行します。現在の試行を終えた後、最後に受理した材料定数・履歴と次の刻みを保存して終了します。新しいフォルダへ`prepare-mc-epoch-recovery`で複製し、複製先の固定実行ファイルで再開します。元ジョブの再実行は拒否します。

`status.json`の`running`は最終報告です。プロセスの生存を保証する値ではありません。異常終了後も`verify-mc-epoch-job`で再開点を検査し、別ジョブへ複製できます。

## 材料経路

仕様は元Groundの材料テーブルと同じ順序で、全材料の目標パラメータを指定します。

```json
{
  "schema": "soil-fem-mc-epoch-path/1",
  "target_parameters": [
    {
      "young_kpa": 30000.0,
      "poisson": 0.25,
      "cohesion_kpa": 20.0,
      "friction_deg": 0.0,
      "dilation_deg": 0.0,
      "tension_cutoff_kpa": null,
      "flow": "non_associated"
    }
  ],
  "initial_increment": 0.0625,
  "minimum_increment": 0.0009765625,
  "maximum_attempts": 128
}
```

E・ν・流れ則の系統・Rankine面の有無を元材料から保持します。cと引張上限を線形、tan(φ)とtan(ψ)を線形で補間し、λ=0を元材料、λ=1を目標材料とします。λはこの補間座標です。

cをc/F、tan(φ)をtan(φ)/Fへ変更すると、この経路の有効低減係数は`1/(1-λ+λ/F)`です。Fは入力した試験条件であり、計算から求めた安全率ではありません。各Fのジョブは同じ未低減施工Runから独立に準備してください。

最小刻み・試行数上限への到達は数値未完です。不収束や材料更新の失敗から崩壊・安全率上端を推定しません。この材料継続ジョブの`stability_classification`は引き続き`not_evaluated`です。要求率との照合は、同じ未低減Runを使う別経路の[epoch SRM](EPOCH_SRM_SEARCH.ja.md)と[工程ゲート](SRM_STAGE_GATE.ja.md)です。一般斜面の安全率やNative UIの完成ではありません。

仕様を変更する場合は、元施工Runから新規ジョブを準備します。封印した保存点の仕様や実行ファイルIDを編集して再利用する操作は拒否します。

## 保存・結果

入力の元Run・メッシュ・仕様と同じ実行ファイルを保持してください。保存点は元Run全体・確定起点・元入力・実行ファイルと、現在の材料定数・履歴・λ・次の刻み・試行列をハッシュで結び付けます。再開時に現在材料の降伏、全体釣合い、変位拘束、固定荷重と施工情報も検査します。

2つの保存ファイルへ交互に書き込み、全書込みを完了してから参照先を更新します。書込み失敗では前の完全な保存点を保持します。強制終了中の未受理のNewton試行は保存対象になりません。

結果の現在変位・応力は、`report.json`の`checkpoint.payload.material_checkpoint.payload.committed`にあります。元Runの`integration_points`との積分点対応を保持します。

材料継続による追加変位は、現在の`displacement_m`から元Runの`last_committed.displacement_m`を引きます。zは上向きなので、追加δuzが負なら沈下、正なら隆起です。表示RESETで構成則用の変位を零へ変更しません。

現在材料は活動要素へ適用します。除去済み要素の履歴は元Runのまま保持し、この材料継続中に新たなBIRTH工程を追加しません。

専用19試験、実プロセスの異常終了後の塑性再開、実CADの施工終点と塑性履歴再開を検証しています。分離環境の全300 Rust/171 Pythonは合格しました。本体実行版と新しいWindows CIの検証記録も確認してください。一般3D塑性の独立解、2D比較、構造/水理との継続と安定判定も未完です。


## 3D結果の出力

```text
python tools/export_mc_epoch.py JOB NEW_OUTPUT
```

最後に受理した保存点を、ジョブ内の固定実行ファイルの`verify-mc-epoch-job`で再検査します。実行中に保存点が進んだ場合は最大3回読み直し、異なる世代を混ぜません。元ジョブを変更せず、その外側の新規フォルダへ出力します。

`origin.vtu`と`current.vtu`をParaView等で開けます。`epoch_displacement_increment_m`は未低減施工終点からの追加変位、`total_displacement_from_seed_m`は累積変位、`displacement_m`は表示基準を差し引いた値です。z正が隆起、負が沈下です。表示RESETでこれらを混同しません。VTU座標自体は未変形のm単位です。

`*-points.vtu`は活動積分点の生の有効応力・塑性ひずみ・間隙比・現在材料定数です。セル平均は`volume_average_`、元節点/要素/積分点IDを保持し、応力を節点へ平滑化しません。引張上限が無い場合はNaNと`has_tension_cutoff=0`で表し、零強度へ置き換えません。除去済み要素はVTUから除き、凍結された全履歴を`accepted-snapshot.json`へそのまま保持します。

途中保存や数値未完でも、検証された最後の受理状態を出力できます。`export-verification.json`に`complete`、到達λ、停止理由、検証世代と元ファイルのハッシュを記録します。数値未完を崩壊判定に読み替えません。2端点の出力で、保存していない中間履歴は再構成しません。

この専用形式は通常の施工Runへ変換しません。Native SolidCraftの専用読込・表示接続は未実装です。出力専用9試験と実CAD6ケース、改ざん拒否5ケースを検証し、Python全180件が合格しました。Rustソースと本体実行ファイルは変更していません。
