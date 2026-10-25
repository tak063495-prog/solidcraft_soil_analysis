# 外部解析ジョブの起動・状態・停止

2026-10-24。拡張版SolidCraftのAnalysisタブから、外部FEMワーカーの起動・監視・停止とCAD変更後の結果失効を実装し、実GUIと実MCC解析で限定検証した。**協調中断・保存点からの画面内再開・3D結果描画は未完で、P5全体の完了ではない。** Python/Gmsh/PowerShellは外部解析側に置き、SolidCraftアプリへ組み込まない。

## 外部CLI

2026-10-25追記。以下の既存nativeワーカー経路と別に、Rustの施工力学CLIへライブ反復通知・受理増分での協調中断・新しいジョブへの再開を追加した。**nativeのStopボタンやONE_WAY/Richards専用コントローラへの接続は未完**である。

新CLIは開発ビルドで実行する。prepareは入力3ファイルと実行バイナリを新規フォルダーへ固定する。run/resumeにはそのフォルダーの固定バイナリを使用する。

```powershell
$solver = '.\target\control-development\release\soil-fem-3d.exe'
& $solver prepare-stages-job MESH.json GROUND.json STAGES.json NEW_JOB
& '.\NEW_JOB\solver\soil-fem-3d.exe' run-stages-job NEW_JOB

# 別の端末から、実際の最終報告を読む／同じジョブへ中断を要求する。
& '.\NEW_JOB\solver\soil-fem-3d.exe' stages-job-status NEW_JOB
& '.\NEW_JOB\solver\soil-fem-3d.exe' stop-stages-job NEW_JOB

# checkpointedを確認した後、別の新規フォルダーへ再開する。
& '.\NEW_JOB\solver\soil-fem-3d.exe' prepare-stages-resume NEW_JOB NEW_RESTART
& '.\NEW_RESTART\solver\soil-fem-3d.exe' run-stages-job NEW_RESTART
```

`progress.jsonl`はジョブID・連番と、実際のNewton反復、物理力残差、混合形式の体積/圧力残差、材料分割数、線形求解診断、受理/棄却を保持する。接触の補助予測子の線形残差はその旨を区別する。状態はprepared / running / checkpointed / succeeded / failedで、status.jsonは原子的に更新する。状態/PIDの記録だけから将来の生存を認定せず、ホストは自分の子プロセスハンドルも観測する。

停止要求は同じジョブIDだけを有効にし、反復途中や棄却後に届いても、新しく受理した増分まで進んでから全履歴・解放荷重・施工継続情報をcheckpoint.jsonへ保存する。初期平衡の途中/終端も対象とする。元の入力/実行版を変えた再開、既存フォルダーへの上書き、同じジョブの二重実行を拒否する。新しい再開ジョブへ古い停止要求を引き継がない。終了値は完了0、失敗1、保存して中断2。checkpointedには完了結果result.jsonを公開しない。数値的に受理できない場合はfailedとなり、中断要求だけを根拠に保存点を捏造しない。強制終了・電源断で新しい保存点が得られるという保証はない。

通知先の失敗はruntime controlエラーとして返し、物理的な刻みの切戻しや破壊判定へ変換しない。観測なしの既存APIと結果の物理フィールドを保つ。新4制御試験・5CLI試験に加え、実SolidCraft由来の斜面で**生存中の施工Newtonを観測してから**停止を要求し、塑性履歴を持つ受理境界で保存、別プロセスで載除荷・局所掘削まで再開した。旧固定版と実行時間/バイナリSHA以外の全結果が一致した（`runs/stage-job-cad-verification-20261025-v2/verification.json`）。初回v1は塑性化前に停止して試験条件を満たさず、その保存点とログ・当時の検証ソースを保持した。同じ開発版の全238 release試験・clippy・書式も合格した。外部CLIの限定検証であり、nativeの画面内中断/再開を含むP5全体は完了扱いにしない。

```powershell
# ワーカーは前景で動く。ホスト側は子プロセスとして起動し、ハンドルを保持する。
.venv\Scripts\python.exe tools\native_job.py run D:\Work\MY-CAPTURE D:\Work\MY-NEW-JOB --mesh-size-m 4

# 別のターミナル/ホストから最後の状態報告を読む、または停止要求を書く。
.venv\Scripts\python.exe tools\native_job.py status D:\Work\MY-NEW-JOB
.venv\Scripts\python.exe tools\native_job.py cancel D:\Work\MY-NEW-JOB
```

`run`には`--solver`、`--order 1|2`、`--quadrature 3|4|5`、JOINTの`--split`を指定できる。ホストは現在の`analysis.soil.snapshot`を確認し、そのSHAを`--expected-input-sha256`へ必ず渡す構成にする。CLIの版不一致は出力フォルダーの作成やプロセス起動前に拒否する。

ジョブ先は新規フォルダーのみ。元のcapture内には作成しない。終了コードは成功0、失敗1、停止2。PowerShellの外側のシェルが非ゼロを1として返す場合もあるため、ワーカー自身のログの終了値と`status.json`を区別する。

## 固定する入力

`capture/`へ全ファイルと完了マーカーをコピーし、コピー前後のハッシュとスナップショットを検査する。`solver/`へ実行バイナリを固定する。JOINT指定もコピーする。CADや元captureが後で変わっても、実行中の物性・形状・工程を差し替えない。

`job.json`には入力SHA、出力元SolidCraft、ソルバー、メッシュ指定と使用するFEMスクリプトのSHAを保存する。実行中にFEMスクリプトが変わった場合は、その結果を成功として公開しない。CADメタデータから実行プログラムを選ばず、ワーカー自身のプロジェクトの既定スクリプトを引数配列で起動する。

内部の`Invoke-SolidCraftCapture.ps1`は`workflow/`に従来どおりのメッシュ・解析・独立検算・VTU・表示データを作る。ワーカーは完了した表示データの入力SHAと実行版、照合結果を確認した後だけ、成功結果への相対パスを公開する。

## 状態報告と停止の意味

`status.json`は一時ファイルからの置換で更新する。状態はpreparing / running / cancelling / cancelled / succeeded / failed。工程は入力固定・形状取込・メッシュ・求解・検算・表示データ準備を区別し、経過時間、ログ末尾、終了値を含める。全体進捗率やNewton反復数を推定して表示しない。

これは最後に報告した状態であり、ファイルがあるだけではプロセスの生存を保証しない。追加したNativeホストは、起動時に得た子プロセスのハンドルを保持し、終了と状態更新を併せて観測する。CAD/属性変更後の結果失効もホスト側で表示する。

`cancel`はジョブIDに対応する要求を書くだけで、保存されたPIDを直接終了しない。動作中のワーカーが自身で起動して保持している子プロセスと、その配下を停止する。Windowsの実プロセスで検証した。既に終了したジョブへの要求は状態を変えない。

**現在の停止はプロセス終了であり、受理増分境界へ到達して新しい再開点を作る協調的な中断ではない。** 未完出力は保持し、停止を成功や解析完了とは扱わない。Rustソルバーの保存再開機能は別に検証済みだが、この停止操作から新しい再開点が必ず得られるとは約束しない。協調中断・ライブ反復診断・保存点を指定したホスト内再開は後続要件に残す。

## 外部ワーカー初回受入の記録

- 単体5試験：異なる入力版の起動前拒否、既存出力と固定入力の保全、停止要求のID・重複・終了済み処理、実際の子/孫プロセスの停止。
- 新しいSolidCraft修正版のcaptureを使ったMCCジョブが、メッシュ・独立検算・VTU/表示データ照合まで正常完了した。最終確定状態と全フレームは、保存した旧ネイティブ版の解析結果と完全一致した。
- 同じ実機出力を細かいメッシュで実行し、求解中に停止要求を送った。ワーカーは停止終了値2とcancelledを報告し、受入済み結果を公開しなかった。停止後、対象ジョブ配下のRustソルバープロセスが残らないことも確認した。

初回の実行記録は`runs/solidcraft-native-ui-fix-20261024-v2/job-success/`、`job-cancel/`、試験ログは`runs/native-job-tests-20261024-v3.log`。当時のネイティブUIとワーカーの対応は`runs/solidcraft-ui-worker-foundation-20261024-v2/`へ保存した。既存の[P6未達要件](P6_STATUS.ja.md)は保持する。

## ネイティブ起動・監視・停止の限定受入

`runs/solidcraft-analysis-source-20261024-v7/`の38ファイル差分に、`analysis.soil.job.start / status / cancel`とSoil AnalysisのAnalysisタブを追加した。Windows実行版のSHA256は`2619f8964813b89e5a89698d2be16799976276961fafafcd08f148d1ff754086`。`runs/solidcraft-soil-jobs-native-20261024-v2/native-binary/`へ実行版を固定した。

専用の起動スクリプトは、既存のポータブル版を使う`Start-SolidCraft.ps1`と分けている。

```powershell
.\scripts\Start-SolidCraftSoil.ps1 -Port 7894
```

このスクリプトは拡張版の開発ビルドを起動し、その子アプリだけに`SOLIDCRAFT_SOIL_FEM_ROOT`を渡す。ルートには`.venv/Scripts/python.exe`、`tools/native_job.py`、`target/release/soil-fem-3d.exe`が必要。アプリは起動時の明示設定からホストを登録し、CADメタデータから実行先を登録しない。計算用のscratchセッションには実行ホストを渡さない。

`job.start`は新規出力先、現在の入力SHA、メッシュ寸法m、要素次数1/2、積分次数3/4/5を受け取る。出力先の`capture/`へ実STEPを出し、ワーカーは`job/`へ別の固定入力とソルバーを作る。同じセッションでは同時に1ジョブだけを起動する。途中の失敗でも診断用のファイルを保持し、既存フォルダーへ上書きしない。

状態報告は起動した`Child`ハンドルの終了と、入力SHA・起動トークン・ワーカージョブIDを照合する。WindowsのvenvランチャーとPython本体はPIDが異なるため、報告PIDとホストの子PIDの一致を前提にしない。トークンはアプリ全体の連番も含め、同じ時計値で別セッションが同時起動しても区別する。完了には実プロセスの終了コード0、ワーカーの成功報告、表示マニフェストと独立VTU照合のSHA・入力版・実行版の一致を必要とする。CAD変更後も固定入力の計算は継続するが、現モデルへの適用を失効させる。Undoで入力が厳密に戻れば適用可能へ戻る。

ネイティブホストは`--watch-stdin`と起動トークンを渡し、専用パイプを保持する。停止要求、所有セッションの破棄、アプリ終了でパイプが閉じると、ワーカーが自身のプロセス木を停止する。PowerShellと解析にはこのパイプを継承させず標準入力を閉じる。プロセス終了後のPIDを再利用して停止しない。停止は引き続き、新しい保存点を作る協調中断ではない。

最新ワーカーでは、実CAD入力の**求解中**にホストパイプを閉じ、終了コード2・cancelled・受入結果なしと、対象ジョブの実行ファイルが残らないことを確認した。`runs/solidcraft-soil-jobs-host-20261024-v3/`へ記録している。接続を開いたままの正常終了も実解析と3回の終了回帰試験で確認した。

実機では2層K0/MCCの局所掘削をAPIとGUIのStartボタンからそれぞれ起動した。両方で旧ネイティブ版と確定状態・全フレームが完全一致し、GUI出力の全STEPバイトも一致した。[GUI正常完了](../runs/solidcraft-soil-jobs-native-20261024-v2/gui-success.png)。別の細かいメッシュをGUIから起動し、ネイティブの所有ハンドルが求解中であることを確認してからStopをクリックした。終了コード2・結果公開なし・残留プロセスなしに合格した。[GUI停止後](../runs/solidcraft-soil-jobs-native-20261024-v2/gui-cancelled.png)。稼働画面は最大0.5秒前の工程表示を含み、`gui-running.png`はメッシュ工程を表示している。停止直前の求解状態は`observed-live-before-ui-stop.json`と実画面操作で別に確認した。

CAD改名中と解析完了後の失効、Undoでの適用回復、異なる入力版・同時実行・既存出力・別ジョブ停止の拒否、終了済み停止の冪等性、不正MCC定数（λ=κ）の失敗処理も確認した。[失効表示](../runs/solidcraft-soil-jobs-native-20261024-v2/stale-result.png)。公開例題の入力を確認してから実際の`app.quit`を求解中に行い、ワーカーが接続切断で停止し、対象プロセスが残らないことも確認した。観測スクリプトの初期報告なし処理を直した際も、同じ実行を観測し直し、ジョブを再投入していない。

Rustの所有権/版/パス/同時IDの7試験、外部ワーカー7試験、既存UI3試験、警告禁止clippy、書式、WASM型チェックが合格した。WASMの既存3ファイルの警告は保持する。最初に発見したパイプ継承待機、Windows特殊パスによるスクリプト拒否、Pythonの終了時バッファーロックは修正し、失敗ログと当時の版を残した。

起動直後、最初の状態ファイルができるまで、一時的な「パスが見つからない」診断を表示する場合がある。状態はrunningのままで次の報告で消えるが、起動待ち表示の整理は残る。過去ジョブの再読込、ライブ反復診断、協調中断/チェックポイント、結果描画、部材/水理を含む全入力方式の画面統合も未完である。
