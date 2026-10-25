# 外部解析ジョブの起動・状態・停止

拡張版SolidCraftのAnalysisタブから、外部FEMワーカーの起動・監視・停止とCAD変更後の結果失効を実装し、実GUIと実MCC解析で限定検証した。Python/Gmsh/PowerShellは外部解析側に置き、SolidCraftアプリへ組み込まない。P5全体は未完了である。

## 2026-10-25：協調ワーカーと再開の接続

`tools/native_job.py run --cooperative` は形状とメッシュを準備した後、Rustの施工ジョブCLIを使用する。停止要求とホストパイプのEOFを、求解中は受理増分での保存要求へ変換する。結果書出し中の停止では完了した生結果を保存し、新しい出力先へ再開して書出しと監査をやり直す。準備中の停止には再開点を約束しない。従来のフラグなし実行は旧プロセス停止方式である。

```powershell
$solver = '.\runs\projected-mcc-development-20261025-v1\solver\soil-fem-3d.exe'
.venv\Scripts\python.exe tools\native_job.py run CAPTURE NEW_JOB --solver $solver --cooperative --mesh-size-m 4
.venv\Scripts\python.exe tools\native_job.py cancel NEW_JOB
# 終了値3とcheckpointedを確認してから、新しい出力先へ再開する。
.venv\Scripts\python.exe tools\native_job.py resume NEW_JOB NEW_RESTART --solver $solver
```

ワーカーの終了値は成功0、失敗1、保存なしの取消2、保存済みの中断3。内部Rust CLIの中断値2と区別する。保存点の入力SHA・ソルバー・処理コード版・データ指紋を確認し、古いジョブから実行プログラムを選ばない。再開時には再メッシュせず、保存時の形状・物性・工程と監査用の来歴をコピーする。異なるCAD入力、異なるインストール済みソルバーや処理版、変更した保存点での再開を拒否する。

実行開始時に外部処理のソースを `runtime/` へ固定し、その版から後続の処理を読む。元作業フォルダーの編集と実行中の処理を混在させない。再開には現在のインストールが保存時と同じ版であることを必要とする。Windowsの一時的な読取競合では、完成済みの同じ一時ファイルを使って置換を最大5回だけ再試行する。拒否が続けば元の状態ファイルを保持して失敗を報告する。

`runs/native-cooperative-verification-20261025-v4/verification.json` に実CADの二層MCC・局所掘削を記録した。連続実行に対して、求解中の停止要求、求解中のホスト切断、結果書出し中の停止の3経路を別プロセスで再開し、`elapsed_seconds` とバイナリSHA以外の全結果が厳密に一致した。全経路の実行版SHAは別に照合している。この例では塑性化を必須にしておらず、塑性履歴の施工CLI検証は後述の別の記録である。単体12試験では入力/保存点の改変拒否、ソース固定、所有プロセスの停止、読取競合も検査した。v2の進捗観測取り逃しとv3のWindows置換失敗は保持している。

SolidCraft側には `job.start` の `cooperative`、`analysis.soil.job.resume`、保存点の指紋検査、進捗表示と再開操作を追加した。旧ホストの7回帰に保存・再開の2試験を加えた9件が合格した。実際の所有プロセスが終了するまで保存状態を認定しないこと、CAD変更による再開の失効、改変した保存点と入力版不一致を出力先の作成前に拒否すること、コマンド経由のプログラム指定の拒否を確認した。記録は `runs/native-cooperative-ui-20261025-v1/engine-verification.json`。ソルバーの選択は起動時の `SOLIDCRAFT_SOIL_FEM_SOLVER` または `Start-SolidCraftSoil.ps1 -SolverPath` だけで行い、CAD入力や保存ジョブの実行パス指定は受け付けない。P1体積MCCの新しい出力形式をビューアへ接続する作業も別途残る。

## 2026-10-26：実画面での停止・再開と別起動からの再開

新しいネイティブ実行版 `5f84d89bcf624c673988146e6e46d21998ae659f4b305b6429b52e7b29da3cdb` をビルドし、実ウィンドウからStart、求解中のStop、保存点からのResume、結果ビューアを操作した。2層K0/MCC・550節点・7776積分点の局所掘削で、進捗0.828125の受理済み増分を保存した。CAD改名でResumeが無効になり、登録コマンドのUndoで入力が厳密に戻ると有効になった。元のモデルを保全した検証用コピーだけを使用している。

同じ起動中の再開と、アプリを終了して保存済みモデルを別起動で開いてからの再開が正常完了した。両方とも同じ入力で128増分を連続実行した結果に対し、全節点・全積分点履歴を含む非実行時フィールドが厳密に一致し、実行版SHAも別に一致確認した。別起動では新規解析用メッシュ欄を既定2 mのままにしたが、再開は保存時の4 mメッシュと入力をバイト単位で保持し、再メッシュしていない。画面のメッシュ指定は新規Start用である。

再開後の実画面で掘削形状、鉛直変位コンター、量とmm単位、完了・監査済みの表示を確認した。[保存状態](../runs/native-cooperative-gui-20261025-v1/saved.png)／[再開後の3D表示](../runs/native-cooperative-gui-20261025-v1/resumed-view.png)。既存UIの編集・CAD変更・入力読込の3回帰試験も合格した。記録は `runs/native-cooperative-gui-verification-20261026-v1/verification.json`。これは従来の変位型MCC形式の連携検証であり、塑性化を必須にした独立材料則検証、P1体積の新形式、全P5/P6の受入を代替しない。

別試験では実際のアプリを新規起動し、登録制御APIから協調Start→求解中の `app.quit` を実行した。終了直前に生存中のワーカーとソルバーへ読取専用のWin32ハンドルを開き、同じハンドルで終了値3と2を確認した。アプリは正常終了0で、進捗0.1328125の受理済み増分を保存した。別起動の実アプリから登録Resumeを実行して正常完了し、128増分の連続解析に対して全非実行時フィールドと元メッシュ・入力が厳密に一致した。アプリ終了後に過去PIDへ停止要求を出す方法は使用していない。記録は `runs/native-cooperative-app-exit-20261026-v2/verification.json`。これは実プロセスと登録コマンドの試験であり、GUIクリック試験とは区別する。v1は試験用出力先を先に作ったためStartが拒否された記録と、観測用Pythonの文字コード問題を保持し、当該アプリを終了して新しいv2で検証した。

検証では固定した開発CLI `09eed93197f8aedcbae79b9d404ccab6becb8b5a0d8b95135cef623a9578210c` をホストのソルバー設定へ指定し、既定の旧受入版は置換していない。総合 `cargo xtask ci` は継続中である。初回はWindowsの `python3` エイリアスでアイコン検査が起動できなかったため、Windowsだけ `python` を選ぶ起動処理を修正した。次の実行では省略取得されていた上流の `skills/` 文書がMCP試験のコンパイルに必要と判明し、同じ固定コミットから `git sparse-checkout add skills` で展開した。試験ソースは変更していない。3回目は書式・資産・アイコン17試験・ライセンス・カタログ整列・依存方向・UIコンテナー・全ターゲットの警告禁止clippyを通過し、全体試験のビルドとWASM等を継続している。ログは `runs/native-cooperative-ci-20261026-v3.log`。未完の品質検査を合格へ変更しない。

以下の2026-10-24の記録と「新しい保存点を作らない停止」の説明は、旧プロファイルと当時のネイティブ実行版に対するものである。

## 外部CLI

2026-10-25時点のCLI追加記録。以下は外部施工力学CLIのライブ反復通知・受理増分での協調中断・新しいジョブへの再開の仕様である。nativeの協調Stop/Resumeとアプリ終了は上記の後続試験で接続・限定検証した。ONE_WAY/Richards専用コントローラへの接続は未完である。

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

## 協調ジョブと新方式に関する後続更新

上記は初回の強制停止型の受入記録である。後続の協調停止・保存再開と新P1体積MCCは [新方式の検証記録](SOLIDCRAFT_PROJECTED_FIELDS.ja.md) を参照する。新方式の実GUIから掘削応力解放の途中で保存・再開し、連続計算の全非時間項目との一致を確認した。自動完了読込のviewファミリー判定修正は別途検証中であり、手動読込の合格で代用しない。

Windowsでソルバーが `status.json` を置換する際の一時的な読込拒否により、ワーカーが生存中の所有ソルバーを停止する問題を再現した。`native_stage_profile.read_solver_status` は同じ256 KiB上限の読込だけを最大5回、10 ms間隔で再試行する。無効JSON、非有限値、存在しないファイルは再試行せず、古い状態を成功として返さない。恒久拒否は失敗として所有する子だけを停止し、出力パイプを閉じる。新3回帰と既存12件の計15試験が合格した。記録は `runs/native-status-read-repair-20261026-v3/integrated-native-job-tests.log`。実再開検証時は同じ修正を隔離配置しており、新規マシンへの配布受入ではない。
