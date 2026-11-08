# 5 kPa低プレロードNC自由土塊：全上部撤去の途中記録

2026-11-08。起点はPR #13の`0136c27bd94047351f6b4cd39bc907664aa31cb7`。
同じ`codex/plastic-spatial-convergence`上の追加調査であり、**この条件の受入は未達**。
P3b-04、構造・水・JOINT全体、P6全体を完了にしない。

## 既存証拠と今回の入力

指定の`FREE_CONTACT`、`JOINTS_VERIFICATION`、`JOINT_SUPPORT`、`STRUCTURE_WATER_JOINT_REMAINDERS`、`P6_STATUS`の関連箇所を棚卸しした。

|条件|既存記録と今回の扱い|
|---|---|
|NC種応力p=pc=100 kPa、追加プレロード100 kPa|自由土塊7条件のうちの限定受入。pcが171.52〜193.62 kPaに硬化した後の除荷・撤去。既存20地点再開を保持する。独立全体変位は弾性MC地盤の範囲で、MCC塑性3D全場の参照ではない|
|NC種応力p=pc=100 kPa、追加プレロード5 kPa|旧長時間試行は中断・未検証。元の失敗入力の生ファイルは今回の既存`runs/`・`interim/`から確認できなかった。旧試行のハッシュ一致を主張しない。今回の新規入力と旧試行を同一実行と扱わない|
|今回の新規5 kPa入力|既存`tools/make_free_contact.py --mcc --pc 100 --preload-kpa 5 --force 1 .5 0`の生成値を使用。54節点・12 Tet10・324材料GP・32接触GP。下部底面のXYZだけ拘束し、上部に拘束を加えない|

初期外圧100 kPaと、工程の追加プレロード5 kPaは別の荷重である。
種材料・幾何・荷重・摩擦・許容値・最小刻み・材料分割上限・生産CLIの数値経路は変更していない。
ジェネレーターと入力依存元も起点と同一である。
今回のバイトSHA256はmesh `73360b7e5718089aec5391f99d13c3b1c42579440589430e3618cfa69e225fd2`、
ground `4b055d60fb99d3e750b086079485f93dd92eeafed0926f29f922d74c4be2a701`、
stages `a69b6e08d93775d9c9bc67503ace9fd9a17627b1e12145744e79167a5483ead2`。
生成済み旧入力がユーザーPCにある場合は、求解より先にこれらのバイトハッシュとJSON条件を照合する。
改行・シリアライズ差がある場合は両方のハッシュを保存し、入力全体の値の一致も確認する。条件差を丸め差として処理しない。

## 今回新たに確認した事実

既存の定期保存CLIを小ケースだけに適用した。
手動の先行観測は約2分で停止要求後、所有するプロセスを終了し、最後の受理済み保存点を保持した。
自動ドライバの別実行は60秒＋終了猶予15秒で停止し、最後の受理点は撤去進行度0.75。
その別プロセス再開は30秒＋5秒で終了し、新たな受理点は得られなかった。再開入力の封印済み0.75保存点は保持されている。
終了後の`status.json`が`running`でも、そのファイルを稼働の証拠とはしない。終了コード・所有子プロセスの終了・保存点を別記する。

先行観測の最高受理点は撤去進行度0.875（λは保持荷重の解放割合）。
受理18増分、完了フレームは初期平衡・載荷・除荷の3枚で、撤去完了フレームはない。
初期平衡・載荷・撤去で意図した収束後棄却を含み、λ=1で`global increment exceeds MCC partition budget`の自然棄却も記録された。
その後のλ=1試行も時間内に受理されなかった。これは**材料更新・数値求解の未完了**であり、崩壊の証明ではない。

|0.875保存点での観測|実測値・検査|
|---|---|
|材料履歴|324点に塑性ひずみ履歴。全324点のpc範囲101.82526050681295〜121.04129818142636 kPa（撤去済み点の凍結値を含む）|
|切戻し|4棄却の確定状態不変。棄却試行に塑性履歴・接触履歴の実変化がある|
|残存域λ=0|接触力の消失を撤去土内力と分けて別再組立した保持荷重と保存監査を、既存のrtol=1e−9/atol=1e−7 kNで比較。自由度の力相対残差3.8035929133263025e−13|
|撤去途中の残存域|保存応力からの別力再組立による相対残差2.442028274979164e−11。既存本体基準1e−8未満|
|凍結|撤去済み上部の全材料GPと全32接触GPの履歴が除荷完了時と完全一致|
|全解放|未到達。保持荷重の残りをゼロと扱わない|

この再組立は**Rustの保存応力・履歴を検査するもの**であり、入力だけから生成した独立塑性3D変位解ではない。
小さい別実行では受理9増分（載荷途中）から新しいプロセスで10増分目へ再開し、初めから10増分まで直接計算した保存点と比較した。
確定状態・全フレーム・地圧初期化/施工の全試行・監査・表示基準が完全一致した。再帰的除外は線形求解の実測時間`elapsed_seconds`だけ。
これは載荷途中1地点の再開であり、全工程途中再開の受入へ読み替えない。

新しい診断ドライバ`tools/probe_low_preload_nc.py`は入力全体を既存5 kPa生成値と照合し、100 kPaケース・入力差を拒否する。
停止要求後も受理境界に戻らなければ、猶予後に自分が起動した子プロセスだけを終了する。
終了済みジョブの最新の封印済み保存点を使用し、再開で新たな受理がなければ元の入力保存点を報告する。
このドライバは受入器ではなく、求解が終了しても独立塑性3D全場・刻み感度・全工程再開を認定できないため`passed: false`を維持する。
観測処理の正常終了コード0は、数値受入の合格を意味しない。

## 未検証項目と受入に必要な証拠

全上部撤去完了と保持荷重ゼロ、入力だけから生成する独立MCC塑性3D全場、両実装の刻み感度、初期平衡・載荷・除荷・撤去それぞれの途中再開と連続経路一致が未完。
今回の途中のλ=0と凍結確認を、全解放・独立塑性変位・メッシュ収束の代用にしない。
既存の弾性MC/100 kPa自由接触参照器を5 kPa MCC塑性の参照器として実行しない。
新しい構成則・亀裂形状・接触BIRTH仕様は作っていない。

受入に進む前に、独立塑性参照が元のMCC構成式・摩擦投影・全体釣合いを入力だけから積分することと、対象全節点3成分・全材料GPの応力/塑性ひずみ/pc/e・全接触GPの力/すべり/仕事の適用範囲を確認する。
Rustの解を参照の初期値・硬化履歴・変位に転用しない。
既存MCC/JOINT検証の元の許容値を保持し、今回の途中出力から新しい参照値や通過数値を作らない。
数値完了だけ、単体試験だけ、静的JSON検査だけでは受け入れない。

## ユーザーPCでの実行

新規出力先を使い、指定起点を含むブランチで実行する。Native UI/GPUは不要。

```powershell
git merge-base --is-ancestor 0136c27bd94047351f6b4cd39bc907664aa31cb7 HEAD
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v

# 5 kPa固定入力＋小さい時間予算での観測。probe.jsonは未検証の診断記録。
.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-low-nc --budget-seconds 60 --grace-seconds 15
# 終了済みジョブから同一exe/入力で別プロセス再開。
.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-low-nc-resumed --resume-from runs/NEW-low-nc/job --budget-seconds 60 --grace-seconds 15
```

時間制限なしで元の数値設定を解く場合は**ユーザーPC**で、既存CLIの定期保存を使う。

```powershell
.venv/Scripts/python.exe tools/make_free_contact.py runs/NEW-low-nc-input --mcc --pc 100 --preload-kpa 5 --force 1 .5 0
& target/release/soil-fem-3d.exe prepare-durable-stages-job runs/NEW-low-nc-input/mesh.json runs/NEW-low-nc-input/ground.json runs/NEW-low-nc-input/stages.json 1 runs/NEW-low-nc-long
& runs/NEW-low-nc-long/solver/soil-fem-3d.exe run-stages-job runs/NEW-low-nc-long
# 別端末から協調停止。受理境界まで戻れなければ停止要求だけでは即時終了しない。
& runs/NEW-low-nc-long/solver/soil-fem-3d.exe stop-stages-job runs/NEW-low-nc-long
# checkpointedの終了ならprepare-stages-resume、異常終了ならprepare-stages-recovery。
& runs/NEW-low-nc-long/solver/soil-fem-3d.exe prepare-stages-recovery runs/NEW-low-nc-long runs/NEW-low-nc-recovered
& runs/NEW-low-nc-recovered/solver/soil-fem-3d.exe run-stages-job runs/NEW-low-nc-recovered
```

全経路が完了した場合だけ、実際の受理列から各工程の途中増分数Nを選び、`checkpoint-stages MESH GROUND STAGES N NEW_CHECKPOINT`と`resume-stages MESH GROUND STAGES CHECKPOINT NEW_RESULT`を別プロセスで実行する。
工程端を途中と呼ばず、全確定状態・全フレーム・全試行/監査を実測時間だけ除いて連続計算と比較する。
再開は同一ソルバーのSHA256を要求する。異なるWindows/Linuxバイナリ間の保存点流用はしない。
刻み感度は元の幾何・材料・荷重・許容値を固定した別フォルダで施工/材料/接触刻みを細分し、独立塑性参照側も細分する。粗側の不合格は残す。
これらの未実行の比較を、この手順の記載だけで合格にしない。

## 保存する既存未達・範囲外

起点の[構造・水・JOINT記録](STRUCTURE_WATER_JOINT_REMAINDERS.ja.md)と`interim/2026-11-08-structure-water-joint/`は保持する。
湿潤JOINTの変位相対L2 `2.6845491097390954e−11`は別条件の合格で、今回の塑性変位の参照ではない。
Fs元区間`[1.3250000000000002, 1.3531250000000001]`、6 Tet10の別区間`[1.3320312500000002, 1.3337890625000002]`、Griffiths中央面L2 15.582／15.886／15.593／15.672％の不合格と連続側判定不能、柱Tet4最細6.828％不合格と連続推定1.034％合格を保持する。
既存MCC掘削の空間差、押出し応力・振動モニター・有限幅奥行き応力の判定不能、杭6.687→4.717％、旧降雨5項目未達と元しきい値も保持する。
実断面SRMは未実行。Fs中点・不収束・試行上限・時間切れを破壊や合格へ変換しない。
指定された大規模/数時間級の再求解は行わない。
部分亀裂先端・接触BIRTH・Richards・板の補間整合などは今回進めていない。

生入力・保存点・観測ログ・試験件数・Windows CI・起点保存照合は[今回の記録](../interim/2026-11-08-low-preload-nc/README.ja.md)を参照する。
