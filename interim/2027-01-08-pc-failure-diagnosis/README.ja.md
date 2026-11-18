# PC失敗6件の切り分けと修正（PR #13、全体未達）

**script bug 1件、environment/runner 3件、genuine numerical/acceptance failure 2件。**
実行中のユーザーPCチェーン（PID 18836）には接続・操作せず、PC計算は実行していない。
PR #13へ積む作業で、マージしない。起点は現在先端
`520cecfce646200505d7fefabaa66ad1bda33f59`。指定`cdaf850c0347c8d841decce7ddff608d03ae7933`の
5コミット後であり、限定depth12/blob:none取得から祖先関係を確認した。
2027-01-08は指定保存先ラベル。PC失敗日は2026-11-18 JST。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6はすべて未達のまま。
今回のソフトウェア修正・限定確認は、元の数値受入の全項目を満たしたことを意味しない。
許容値、参照値、cutback/line search打切り、生産Rust、NCコード/設定は変更していない。
既存診断スクリプトもバイト保持し、新しい入口を`pcfail_*`へ移植した。

## 一次資料と以前の結論を引き継ぐ範囲

PCの根拠は依頼本文の[埋め込みログを転記した一次資料](PC-PRIMARY.ja.md)。
失敗時点の6件のSHA・数値・例外・未到達範囲を保持する。PCのrun.json/生ログそのものを取得・
SHA検証したとはしない。11-15の2GB capフルcheckout失敗ZIPは別件で、今回の一次資料には使わない。
クラウド新規求解とPCログの数値を混ぜず、[measurements.json](measurements.json)に分ける。

参照した既存資料とソースのSHAは[prior-record-binding.json](prior-record-binding.json)に固定した。
今回、外部論文本文を新しく取得したとの主張や、既存文献値を変更したとの主張はない。

- [B-02 12-24](../2026-12-24-b02-offset-composite-cause-diagnosis/README.ja.md)と
  [01-02追補](../2027-01-02-b02-offset-composite-cloud-followup/README.ja.md)：12/96の限定比較、
  Q/仕事の空間差、薄板96の残差未達と768の旧失敗を保持。板厚/腕/宣言順/積分を再診断して
  新しい単一原因にしていない。今回の対象はWindowsのプロセス識別失敗。
- [B-01 12-25](../2026-12-25-b01-curved-tri6-cause-diagnosis/README.ja.md)と
  [01-01追補](../2027-01-01-b01-curved-tri6-cloud-followup/README.ja.md)：ガード省略の内部対照は
  真の曲面要素を受け入れるものではない。共変性・弱いdrilling方向・丸め行列のDecimal対照の旧結果を保持。
  今回の屋根n8停止を、以前の150DOF共変性対照や曲面ガードの失敗と同一視しない。
- [A-04 12-23](../2026-12-23-a04-richards-wet-joint-cause-diagnosis/README.ja.md)、
  [11-29の独立SOIL参照](../2026-11-29-richards-wet-joint-a04/soil_reference.py)、
  [12-31追補](../2026-12-31-a04-richards-wet-joint-cloud-followup/README.ja.md)：水頭種/所有写像・時計・
  連成仕様の不足、h×dt×VG積分、宣言順による乾燥Newton停止の旧失敗は別の残件。
  今回の参照line searchを、PDE離散化誤差や未定義JOINT写像へ置き換えない。
- [A-01 12-20](../2026-12-20-a01-joint-combo-cause-diagnosis/README.ja.md)：仕事差・平行面・
  旧線形停止と全未達を保持。今回の証拠欠落はこれらを再計算した結果ではない。
- [C-04 11-22](../2026-11-22-wet-joint-transient-water/README.ja.md)は独立RT0参照で
  既に丸め停止とfloat128 LAPACK非対応を記録し、局所4流束の拡張精度化で元の厳しい残差を維持していた。
  現行[hydraulic_reference.py](../2026-11-22-wet-joint-transient-water/hydraulic_reference.py)の
  `retention/local_solve/Model.step`を直接確認した。拡張精度がどのOSでも同じとは限らない点を今回交差した。
  [01-03追補](../2027-01-03-c04-transient-water-cloud-followup/README.ja.md)で消化済みの16/64 P/R交換を
  再びクラウド残件へ戻さない。一般化の仕様・実装待ちは未達のまま。

## 最終分類と根拠

|項目|主分類・主因|副因と切り分け|今回の到達点|
|---|---|---|---|
|b02-n4|**script bug**：終了済みnativeと別のresumeにPID不一致を要求|Windows PID再利用は正当。全該当子は終了0、26.6s/約125MiB WSでありtimeout/memcapではない。exact比較前のassert停止|新入口で起動記録をUUID化、Run/後続状態のexact検査を維持。旧assertを同PIDで再現する試験あり|
|b01-evidence|**environment/runner**：sparse除外の12-25証拠part不在|欠落を記録する契約に対してp.statが例外停止するscript側の扱いも不足。11試験中の同一欠落1失敗を数値不合格へ転用しない|欠落全partをincompleteで保存、終了2。必要2partのみ取得して全entry SHA検証|
|b01-cli|**environment/runner**：12-03証拠part不在|旧拒否入力の読込前に停止、solver未起動。4partの別依存でb01-evidenceと同一アーカイブではない|事前検証入口、欠落時solverを起動せず記録。選択取得・旧21拒否/4受理の検査を維持|
|b01-roof-n8|**genuine numerical/acceptance failure**：梁/板内部系の丸めCholesky因子による真残差成分が主で、元1e-10を超過|n2/4成功、n8のみcutback枯渇。平面多面体で幾何ガード拒否ではない。固定順別残差も超過し、表示集計だけの問題ではない|同じn2/4/8をクラウド再現。棄却系を保存、同一行列の別求解で丸め感度を限定確認。n8完走未達|
|a04-matrix|**genuine numerical/acceptance failure**：独立RT0参照のline searchが元残差条件を達成できず停止|native成功。副因は参照の精度可搬性：binary64状態評価で丸め停滞を再現。外側RuntimeError・ZIP重複名警告は停止源ではない|元方針のまま状態精度×局所流束精度を交差。PCの元失敗のbit一致は未確認。Windows参照の完走未達|
|a01-evidence|**environment/runner**：12-20と11-23の異なる証拠partが不在|verify/extract-old双方とも求解前。後続Tet4等doneと両立し、A-01数値合格の証拠にはならない|1+5partだけ取得・全entry検証、旧選択展開範囲を維持。欠落はincomplete、終了2|

原因寄与を足して100%に配分しない。peak_ws=0を必要メモリ0と解釈しない。
埋め込みログで示された6停止をtimeout/memcapで説明しない。

## スクリプト・環境の修正

[pcfail_b02.py](pcfail_b02.py)は先端のb02cause_driverを移植した明示実行入口。
`Popen`ごとのinvocation_idとfresh_popenを記録し、OS PIDも観測値として残す。
待機後のPIDの大小/不一致をプロセス同一性の条件にしない。
別起動・終了0に加え、**Run全体とresume後のaccepted状態列は旧exact比較のまま**検査する。
除外は従来の実測時間だけ。状態・応力・受理列の改変は回帰試験で失敗する。
旧driverと旧失敗記録は修正・削除しない。PCでは新入口を使う必要がある。

[pcfail_evidence.py](pcfail_evidence.py)は4群のmanifestに列挙されたpartだけを扱う。
missingは全件列挙して`status=incomplete, passed=false, numerical_acceptance=false`、終了2。
サイズ/SHA不一致はfailed、終了1。verifiedでも数値passedはfalseのまま。
全part・結合ZIP・全entryを検証してから必要範囲だけを展開し、欠落時に空の「展開済み」先を作らない。
固定40桁commitのgit blob選択取得を行い、全履歴取得/全面checkout/既存破損partの上書きをしない。
今回必要な12partは合計207,550,464B、b01展開58,390,662B、a01-oldは旧選択部分だけ。

[pcfail_cli.py](pcfail_cli.py)は旧B-01 controlsを事前検証後に呼ぶ。
[pcfail_run.py](pcfail_run.py)は新規子だけを資源制限して、非ゼロ・incomplete・ログを保存する。
独立項目のchainでは非ゼロを消さず、後続の独立項目も記録する。依存する数値高水準は自動登録しない。
Windows Job経路はLinuxでは実行できず、Windows CI/PCでの確認が残る。

[軽量回帰試験](../../tests/test_pc_failure_diagnosis.py)には旧PID assertの再現、同PIDの二つのPopen、
exact改変拒否、欠落/破損/全SHA、実sparse Gitからのpart選択復元、非ゼロ後の独立続行・timeout、
状態精度の隔離、観測用コピー、実CRLF文書とPOSIX/NULパス比較を含める。
重い数値求解をCIの新試験へ追加しない。rust-integration-1への追加0件。
CI sparseには今回の4partをPython coreだけへ含める明示例外を追加し、他群では除外する。
実workflowの条件を解決し、POSIX相対パス・NUL入出力・LF/CRLFの両規則でGit check-rulesを照合する。
初回全履歴転送や29GB展開を新しいPC手順へ追加しない。
初回に見つかったsparse親ディレクトリ欠落の実装不具合は修正し、初回10成功/1errorのログも残す。
封印追加後の試験では固定順残差関数へlistを渡したfixtureの型誤りも検出した。
実計測側は最初からndarrayであり、試験入力だけを訂正した12成功/1errorの初回ログも保持する。

## b01-roof-n8：測定・離散化・線形求解の区別

[pcfail_beam_build.py](pcfail_beam_build.py)はwork内へソースをコピーし、元の失敗枝でだけ、
環境変数PCFAIL_BEAM_TRACEが明示されたとき最初の棄却系をcreate_newで保存する。
本体ソースは不変。既定オフ、演算・解の採用・許容値・cutbackに変更0。
観測パッチ/元SHA/生成SHA/ビルドログ・生入力・Runを封印する。

クラウドでもn2/4は完走、n8は同じ`load_plate_beam: ... beam interior true residual exceeded tolerance`。
n8の最初の棄却はfactor=.25、1584内部未知数、残差2.948857986787151e-10。
[pcfail_beam_report.py](pcfail_beam_report.py)で物理行/列の順を固定したmath.fsum再組立は
2.918568500898291e-10で、元ゲート1e-10超過が残った。
PCのlast progress .001953125と、この最初の棄却factor .25を同じ保存点として扱わない。

保存行列の対称性差max=5.960464477539063e-8。
下三角を対称に解釈するCholesky対照の固有値端点は約289.463/1.211012599e9、
混在DOFの代数的条件比は約4.18365e6。単位の違う力/モーメントDOFを混ぜた量で、
物理力基準や物理座屈の判定にしない。非対称な全保存配列そのものの厳密固有値とは呼ばない。

同じ丸め済み行列・右辺にSciPy Choleskyを適用した固定順残差は7.981058583503436e-11。
元解への固定順残差補正を観測だけ行うと第1補正で2.2906647606869527e-11。
いずれも保存された**全配列**に対して残差を計算し、下三角だけで評価を偽装しない。
この対照は同じ離散系で線形求解/丸めの寄与がある証拠で、単なる最終ノルム順序差ではない。
元K生成の丸め・固有空間・各拘束/drilling寄与までの完全な分解はしていない。
別解を本体へ入れたり、反復改良を生産ソルバーへ実装したりしていない。

追加の既定オフ・work-only観測で丸めCholesky因子Lと前進代入yも保存した。
元の最初の棄却行列A・右辺b・解xは全成分厳密一致し、同じLからの三角求解再生もxと厳密一致。
`src/beams.rs:955–960`がkqq.clone()をそのままCholeskyへ渡すことと、nalgebra 0.34.2の
前進/後退代入の実装を確認した。対角補正や異なる参照行列を因子誤差へ混ぜていない。

[因子/三角求解の分解](factor-components-observation.json)は、Aの下三角を対称化したSを用い、
`(A−S)x + (Sx−L(Lᵀx)) + (L(Lᵀx)−Ly) + (Ly−b) = Ax−b`の符号付きベクトルを保存する。
[Decimal 50/80桁の独立積和照合](factor-decimal-observation.json)では、同じ丸め係数だけから
固定DOF順に再評価した（新規solver実行0）。ノルムは明示したfloat64変換後の固定順math.fsum。

|Decimal80桁の残差成分 / 右辺L2（分母689.5988761109853）|相対L2|
|---|---:|
|下三角の対称化、(A−S)x|1.421209412023646e-11|
|丸め因子、Sx−L(Lᵀx)|2.79919286743232e-10|
|後退代入を元方程式へ写した成分、L(Lᵀx)−Ly|6.201283292488465e-11|
|前進代入、Ly−b|3.48739655483843e-14|
|全成分、Ax−b|2.931542484607753e-10|

因子/前進/後退代入の各1584成分は50/80桁のfloat64丸め戻しで全成分一致した。
対称化項には12成分・最大1.2924697071141057e-26の差があり、そのまま保存し、
許容値で一致扱いにしない。各成分ノルムのfloat64値は50/80桁で一致した。
丸め因子項だけでも元ゲートを超え、前進代入・対称化の項はその大きさを説明しない。
**この固定された棄却系の主寄与はCholesky因子化側と限定同定できた。**
相関するベクトルなのでノルムは加法的な原因割合ではなく、和が100%となる分解ではない。
追加の通常経路も同じ12試行で未完了。因子化の精度改善を本体へ実装していない。

[切戻しの入力代数](roof-cutback-algebra.json)も保存した。屋根は梁0・土付着0・節点荷重0、
板の定数面荷重/自重0、保存された設置参照も全成分0で、同じ活動配置の線形板系である。
荷重係数.25から.001953125までの棄却8回は全て確定状態不変。
同じ系で右辺だけを二分すると、右辺L2の予測値は最後でも5.387491219617073で元の分母下限1を超える。
したがって一様な荷重縮小だけでは同じ線形系の相対残差を改善できない。
これは保存入力/場からの代数であり、未保存の各試行のnative残差を新しい実測値として作っていない。
切戻し上限を増やして通す修正は行わない。

n2/4→n8ではメッシュと系の条件がともに変わるので、有限要素の不整合や曲面理論の誤りを
単独原因と断定できない。反対に、この丸め対照の成功で屋根の全場収束・真の曲面殻・
文献受入を満たしたとも言えない。**n8の求解は未完了で、修正済みとはしない。**

## a04-matrix：参照精度の交差、PDE誤差との区別

[pcfail_hydraulic.py](pcfail_hydraulic.py)は入力だけから旧SOILモデルを新規に作る。
NumPyそのものをグローバル変更せず、隔離した参照のlongdouble型だけを代理指定する。
元の60Newton、20line-search候補、merit<.01、真線形1e-10、dt=.25、12 Tet4、
VG/Mualem・貯留差・力学/水理境界・積分則を保つ。

|状態・貯留評価の型 / 局所4流束求解の型|完走|停止merit（クラウド新値）|
|---|---|---:|
|host拡張（nmant63）/ host拡張|true|全経路完走、24状態|
|binary64（nmant52）/ binary64|false|0.01900050585170866|
|binary64 / host拡張|false|0.011983007999600172|
|host拡張 / binary64|true|全経路完走、24状態|

24は参照の初期化4＋nativeの20frames相当の合計。PCの`construction complete:20frames`と
矛盾しない。新しい24をPCの保存値へ上書きしない。

binary64停止時の残差L2=1.0889474740275325e-17、最大Newton更新=6.357697777811309e-17。
cut0で18成分、cut1で1成分だけ変化、cut2〜19は状態の変化0。
20候補の最小scaled residualも0.01900050585170866で元条件.01を満たさない。
状態binary64/流束拡張でも残差8.969037500620236e-18、最大更新5.104207261490814e-17で
cut2以降の状態変化0、候補に元条件を満たすものがない。
従ってこの対照の停止は、候補数不足を延長すれば解消するとの証拠ではなく、
状態更新/貯留差の達成精度に関わる丸め停滞である。局所4×4求解だけの精度化では足りない。

[保存された失敗状態のULP観測](precision-floor-observation.json)では、binary64の最悪成分は
残差−2.590768969595179e-18 m³/sに対し、元条件から決まる絶対目標が1.3635263133598091e-18 m³/s。
42未知数中41成分のNewton更新が保存状態の1 ULP未満、24成分が半ULP未満だった。
これは保存データの再観測（追加solver実行0）であり、誤差保証や許容値の新設ではない。
元の厳しい参照判定とbinary64の達成精度が接近していることを、状態が変わらない候補の記録と併せて示す。

PCログの0.011860214796293166と上表は**同じ数値ではない**。
PCのNumPyビルド・finfo・失敗候補場を未取得なので、元PC失敗のbit一致/全算術寄与は未確定。
Windowsでlongdoubleがbinary64相当の場合の説明に強い実験的根拠を与えるが、
未取得のPC状態を確認済みと書かない。再実行手順でfinfoと部分場を保存する。
残差条件を緩和したり、止まった参照を独立正解として比較したりしない。

同じ入力・空間/dt/積分則の精度対照で完走/停止が分かれたため、今回の停止を
h/dt/VG積分のPDE離散化誤差だけで説明できない。
既存12-23/12-31のPDE差、JOINT非対応、所有/物理時計/連成仕様は別途未達。
portableな参照精度の実装・検証は今後のクラウド残件で、Windows matrix完走は未確認。

## 証拠・再現・残件

クラウドの生入力、失敗Run、棄却線形系、4精度対照の全反復・候補・失敗状態、
選択証拠検証、旧PID再現試験と初回失敗、ソフトウェア検証ログを新フォルダ側へ保存する。
重いデータはarchive-index.jsonと分割ZIPの全SHAで封印し、旧封印へ追記しない。
今回の科学証拠は61,109,258B・4part・663entry。全part/結合/entryを再読込してSHA照合した。
B-02のn4は12/96の2経路・40受理状態・5地点の再開が全てexact一致、元限定場ゲートも充足。
観測器は変更のない生産ライブラリから直接rustcでビルドした。driver内scopeはwork-only観測版CLIの
trace無指定で実行し、別途通常の生産CLIでも旧scopeを終了0で確認した。
B-01の通常CLI controlsは21拒否/4受理が全契約一致。A-04の通常CLI対照は20frames/1removalで完走。
これらを屋根n8・Windows独立参照・全受入の完了へ読み替えない。
PC元ログを再取得したとの主張はなく、元PC失敗出力はその場に保持する。

[REPRODUCE.ja.md](REPRODUCE.ja.md)に各cwd、固定SHAのsparse取得、Python3.13 venv/Rust1.99.1、
autocrlf=false、BLAS/MKL/OMP=1、一意出力、予算、見る値、依存条件を記載した。
b02-n8/n16はn4の完走・厳密再開・予算確認後だけ。A04の32→64→128もmatrix成功後だけ。
どちらも旧skippedを自動再実行へ登録しない。768/6144・40GB超・1経路1日超の新規PC計画なし。
必要なら上限を超える一般確認は**規模上限のため確認不可・未達**とする。

|対象|クラウド残件|PC再実行|判断/実装待ち|
|---|---|---|---|
|b02-n4|新入口のWindows可搬性以外の旧B-02残件は01-02を維持|要：新入口でn4。高水準は条件付き|元の接続/設置参照/世代/仕事の仕様を本件で選ばない|
|b01-evidence / cli|必要partの選択取得とSHA検証は本件で確認。Windows Job経路の確認は残る|要：必要partのみ取得して証拠/controls|修正について新しい物理仕様判断なし|
|b01-roof-n8|固定棄却系の因子化/前後代入の分離は消化。K生成・拘束/drillingと条件の関係の一般化は小型で可能。上限外へ移さない|任意：同じ通常CLIの再失敗を保存|生産求解を変える場合は別の実装判断。曲面要素の仕様待ちは別件|
|a04-matrix|portableな状態/貯留/残差精度の参照実装と独立検証、PC元状態との照合|要：matrixとfinfo。32以降は成功後だけ|参照精度実装の選択。水理JOINTの所有/時計/連成仕様は従来どおり別途判断待ち|
|a01-evidence|選択part/全entry/旧subset検証は確認。数値A-01残件は12-20を維持|要：2群の証拠を新入口で確認|新しい物理仕様判断なし|

P6-01へacceptance/measurements/確定validationをそれぞれ固定コミット・Git blob・SHA256・
型付きboolで登録する。台帳の元バイトを接頭部保持し、監査器は台帳SHA定数だけ更新する。
最終件数・全実行/初回失敗・再実行の選択は[validation.json](validation.json)に記録する。
**push後のWindows CIは未確定**。Linux成功をWindows成功へ転用しない。

## Linux全件検証の確定値

fmt・全ターゲットclippy・locked release全ターゲットbuildは終了0。
Rust全ターゲット＋docは467成功、Python tests/toolsは979固有件成功（新規13件・tools2件は内数）。
最終選択の失敗/error・未実行・skip・ignore・filter除外・expected failure・unexpected successは0。
Rustの独立一覧/実行IDと、Pythonの最終発見/選択IDはそれぞれ完全一致した。

[VALIDATION.ja.md](VALIDATION.ja.md)と[validation.json](validation.json)は、最初の部分実行・中断・
欠落blobによる非ゼロ試行・CI固定行数の失敗・修正後再実行を別に数える。
全件の生試行1213 test開始＝984成功/1失敗/228test errorに、test未開始等のfixture error74、
中断1を別記する。これを最終979固有成功へ二重加算しない。
準備/単独回帰の2errorも保持。登録後確認は別記し、全件へ加算しない。

全件ログと追加の因子/ULP観測は12,135,145B・1224entryのfull-suite.zipへ封印し、全entry SHAを照合した。
検証用bindingの1653コード/設定は保存後から終了まで全SHA一致。
P6登録に必要な監査器SHA pinの更新は、登録後の専用試験と最終監査で別に確認する。
**push後のWindows CIは未確定。数値の全受入は未達。**

確定validationは`607c66c89fd9093fd718828a408eba1bf7655777`に固定して同じ台帳へ登録した。
最終監査は出典369件、integrity_ok=true、SOURCE_NOT_INVENTORIED=0。
登録後のP6/CI/新試験92件も全件成功し、467/979へは加算していない。
[最終監査](final-p601-audit.json)、[登録後結果](post-registration.json)、[保全検証](preservation.json)を参照。
旧台帳は元バイト接頭部を保持し、監査器はSHA pinだけの変更。旧不変ファイルはGit objectを保持し、
参照121実ファイルは全SHA一致。29GBを展開して全実ファイルを再hashしたとの主張はない。
PC実行0・PC実行中チェーンへの操作0・生産Rust/NC変更0。push後Windows CIは未確定。

今回の限定確認を全体完了へ変更しない。
