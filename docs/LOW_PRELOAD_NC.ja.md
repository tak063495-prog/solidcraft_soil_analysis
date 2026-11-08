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

## 85556a4からの続行調査：全解放は未達（2026-11-08）

上の前回記録を保持し、指定先端`85556a460179b96ccd95f959cdca7a3b6ff8e156`から同じブランチで続行した。
**最高受理進行度は0.95123291015625、受理259増分。1.0には到達せず、計測JSONは`passed: false`のまま。**
初期平衡・載荷・除荷の3フレームだけが完了している。未収束・切戻し上限・線形残差超過を破壊や合格に変換しない。
P3b/P4の行、構造・水・JOINT全体、P6全体を完了にしない。

### 保存点からの再開と収束制御の比較

元の数値ソースは変更していない。再ビルドしたLinux実行版のSHA256は
`afc5b3f542bd2bd68b9520f5a392780cb23fbf8cdc630f689f71a6ca57dd260e`で、前回0.875保存点の封印と一致した。
元ZIPを展開したコピーへこの同一バイナリを復元し、既存`prepare-stages-recovery`で再開した。
保存点の内容・チェックサム・入力・バイナリ結合を手で書き換えていない。
元の0.125継続刻みから、既存の自然切戻しだけで0.94921875（受理21増分）まで進み、実測521.608264921秒で切戻しを使い切った。
600秒予算より先の数値停止であり、前回の時間切れとも区別する。

封印済み保存点の設定を変更するAPIは既存CLIにないため、設定を変える比較は同じ5 kPa入力を最初から計算した。
新しい診断オプションは撤去の`initial_increment`、既存Newton上限の範囲内の反復数、最小刻みの細分化、既存dense LUへの切替だけを許す。
物理入力全体を元ジェネレーターと比較し、5 kPa、p=pc=100 kPa、水平1/0.5 kN、形状・拘束・構成則・摩擦・力/線形許容値・材料/接触局所刻みを固定する。
dense LU時は既存2000自由度上限に収まる162のDOF予算を指定した。最初に元の疎解法用100000予算を残した2試行は入力検査で拒否され、受理点なしとして保存した。

|解法・撤去初期刻み|Newton上限|最小刻み|最高受理進行度|停止理由|
|---|---:|---:|---:|---|
|元の0.875保存点から、GMRES|60|1/512|0.94921875|Newton上限後、切戻しを使い切り|
|GMRES・1/16、1/32の2条件|12|1/512|両方0.947265625|Newton上限後、切戻しを使い切り|
|GMRES・1/32、1/64の2条件|60|1/512|両方0.94921875|Arnoldi breakdown、真の線形残差が元の1e−10を超過|
|GMRES・1/64|60|1/4096|0.951171875|Newton上限後、切戻しを使い切り|
|dense LU・1/64|60|1/8192|0.951171875|Newton上限後、切戻しを使い切り|
|dense LU・1/128|60|1/16384|0.95123291015625|Newton上限後、切戻しを使い切り|

意図した収束後棄却を保持したので、例えば初期1/128の最初の受理幅は1/256となる。
最高点の最後の棄却残差は`2.8009435309853623e−5`で、元の力基準1e−8を満たさない。
受理した最高点の保存応力からの別再組立残差は`4.114984391472963e−10`、残存域λ=0は`2.4822464644593524e−14`。
撤去済み162材料GPと32接触GPの全履歴は除荷完了時と完全一致し、全324材料GPに塑性履歴がある。
λ=0、共通進行度0.875/0.9375、各経路の最終受理点で全GP・接触の監査を保存した。
全解放後の凍結・保持荷重ゼロは未確認であり、現在の凍結をその代用にしない。

### 独立塑性3D全場の部分チェックと刻み感度

`tools/low_preload_nc_reference.py`は既存の独立MCC連続体速度式・適合Tet10全B行列・中点/RK4イベント積分と、既存JOINTの閉鎖固着枝のばねを独立全体組立で結合する。
入力だけから全節点3成分変位と全324材料GPの応力・塑性ひずみ・pc・e、全32接触GPのジャンプ・力・すべり・弾性エネルギーを計算する。
本体の変位・応力・接線・履歴・保持荷重・受理刻み列を参照の生成に使わない。
各受理参照点で閉鎖・Coulomb固着条件を検査し、開離・滑りが生じればこの限定参照は拒否する。新たな構成則や接触仕様は作らない。
撤去後の非活動材料・接触履歴は独立経路でも保持する。

RK4の128/256分割はそれぞれ撤去`0.9509514227042946` / `0.9509514099202983`まで進み、真の線形残差`4.490642967739375e−10` / `3.582143726547472e−10`が既存1e−10を超えて停止した。
中点128/256も約0.939で既存の降伏面ドリフト制限に達し、投影や許容値緩和で続けていない。
**独立塑性3Dの全撤去経路の完了・受入は未達**。保存応力の釣合い監査と、入力だけから生成したこの部分参照を区別する。

|共通進行度|RK4 128→256の全節点変位相対L2差|本体1/32→1/64の全節点変位相対L2差|本体1/64対RK4-256の全節点変位相対L2差|
|---|---:|---:|---:|
|0.875|1.6645670575507345e−7|0.0011652650119939995|0.007109946364389393|
|0.9375|7.144051321864514e−6|0.0229499175883048|0.013723275397741966|

全材料GPの応力・塑性ひずみ・pc・eと接触場の全点/活動点差も計測JSONへ記録した。
これらは測定値であり、新しい合格数値・参照定数・空間収束の認定ではない。
接触端点仕事・数値散逸・塑性散逸は既存`joint_reference.System.audit`のエネルギー恒等式で検査した。
本体の離散端点仕事と独立連続率参照の仕事を同一視しない。

### 途中再開で今回確認した範囲

初期平衡0.5、載荷0.25、除荷0.5、撤去0.25から別CLIで次の受理点まで再開した4比較は、直接経路と完全一致した。
比較は入力・レイアウト・継続情報・全確定状態・フレーム・全試行/監査を含むRun全体で、再帰的除外は`elapsed_seconds`だけ。
これは各工程の途中比較であり、未到達の撤去完了まで一致したという意味ではない。
細分経路の撤去0.9375からの追加比較も別記録へ保存する。

### ユーザーPCでの続行・再測定

Windows実行版はLinux保存点とバイナリSHAが異なるため、Linux保存点をWindows用に書き換えず、同じ生成条件から計算する。
すべて新しい出力先を使う。各診断の終了コード0は観測処理の完了だけを表し、数値受入ではない。

```powershell
git merge-base --is-ancestor 85556a460179b96ccd95f959cdca7a3b6ff8e156 HEAD
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
.venv/Scripts/python.exe -m unittest discover -s tests -p 'test_*.py' -v
.venv/Scripts/python.exe -m unittest discover -s tools -p 'test_*.py' -v

# 最高到達の数値設定を、同じ5 kPa物理入力で再計算。
.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-nc-dense128 --budget-seconds 600 --grace-seconds 5 --excavation-increment .0078125 --minimum-increment .00006103515625 --dense-lu
# 終了したジョブの最新封印保存点から、同一exe・同一設定で再開。
.venv/Scripts/python.exe tools/probe_low_preload_nc.py --solver target/release/soil-fem-3d.exe --output runs/NEW-nc-dense128-recovered --resume-from runs/NEW-nc-dense128/job --budget-seconds 600 --grace-seconds 5 --excavation-increment .0078125 --minimum-increment .00006103515625 --dense-lu

# 独立参照の2刻み。途中拒否もreference.jsonにpassed=falseで残す。
.venv/Scripts/python.exe tools/low_preload_nc_reference.py --output runs/NEW-nc-reference128 --steps 128 --method rk4 --budget-seconds 600
.venv/Scripts/python.exe tools/low_preload_nc_reference.py --output runs/NEW-nc-reference256 --steps 256 --method rk4 --budget-seconds 600
# 全4工程途中の再開比較。全撤去合格を表すコマンドではない。
.venv/Scripts/python.exe tools/check_low_preload_nc_restarts.py --solver target/release/soil-fem-3d.exe --output runs/NEW-nc-phase-restarts
```

同じ数値設定で切戻しを使い切った場合、実行時間だけを延ばせば進むとは主張しない。
さらに最小刻みを細分する比較は、同一物理条件を固定した別フォルダで最初から実行し、粗側の停止記録を保持する。
設定を変えた入力へ古い封印保存点を付け替えない。力基準・線形基準・MCC分割上限を緩めて合格にしない。
全解放、独立全経路の完了と刻み感度、撤去完了までの再開一致が揃うまで`passed=false`を維持する。

生証拠、設定差、全場計測、今回の全試験件数と先端Windows CIは[続行記録](../interim/2026-11-08-low-preload-nc-followup/README.ja.md)に記載する。
