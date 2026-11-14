# C-03：塑性バルクと湿潤JOINTの原因切り分け

**C-03 は未達、passed=false。** 限定確認は小型離散問題の診断・全履歴比較・実再開に限る。一般バルク塑性との湿潤統合を完了にしない。C-01・C-02・C-04、NC-01、P6、A/B各項も未達のまま。

PR #13、`codex/plastic-spatial-convergence`、起点 `1802c11be0145c706c49d4b698b2e5aed235c509`。マージしない。保存先2026-12-18は依頼された作業系列ラベル、実施2026-11-14（日本時間）。[前回C-03](../2026-11-20-wet-joint-plastic-bulk-integration/README.ja.md)、[C-01原因診断](../2026-12-13-c01-cause-diagnosis/README.ja.md)、[C-01全体経路対照](../2026-12-16-c01-regen-global-control/README.ja.md)、[C-02](../2026-12-17-c02-spatial-cause-diagnosis/README.ja.md)を参照した。

生産 `src/`、旧入力・旧JSON・旧hash/snapshot・旧オプション・旧失敗記録は変えない。新しい `c03_bulk_cause_probe` は既存Observerへコピーを要求する `--diagnose` のみで、既定オフ。旧exampleも変更しない。新Python名は `c03diag_*` とし、既存 import 名へ割り込まない。NCコード・設定は変更しない。todoの既存C-03行は未達のまま保持する。

## 先行事例（実験より先に確認）

確認日: 2026-11-14（環境の日付。フォルダ名は依頼された作業系列名）。以下は本文を確認できた範囲の要約で、原文の引用ではない。

- S. W. Sloan (1987), *Substepping schemes for the numerical integration of elastoplastic stress–strain relations*, IJNME 24, 893–911, DOI 10.1002/nme.1620240505。[著者所属大学の本文](https://www.newcastle.edu.au/__data/assets/pdf_file/0017/22463/12_Substepping-schemes-for-the-numerical-integration-of-elastoplastic-stress-strain-relations.pdf)、pp.893–894: 既知のひずみ増分に対し modified Euler / Runge–Kutta と局所誤差推定で内分割する。局所積分の検証と境界値問題の検証は別。今回: 端点固定の局所再生は局所感度の対照に限定する。
- S. W. Sloan, A. J. Abbo, D. Sheng (2001), *Refined explicit integration of elastoplastic models with automatic error control*, Engineering Computations 18(1/2),121–154。[大学公開本文](https://www.newcastle.edu.au/__data/assets/pdf_file/0006/22578/47_a_Refined-explicit-integration-of-elastoplastic-models-with-automatic-error-control.pdf)、pp.121–123 を確認。既知経路に沿うひずみ内分割の誤差制御、陰的 return mapping / 整合接線 / 頂点の扱いを区別している。今回: 材料と接触の分割を独立に変え、全体経路を別に比較する。2002 erratum の[存在とPDF](https://www.newcastle.edu.au/__data/assets/pdf_file/0007/22579/47_b_Erratum-Refined-explicit-integration-of-elastoplastic-soil-models-with-automatic-error-control.pdf)も確認したが、文字抽出できず訂正式の内容は未確認。論文のアルゴリズムを新規実装する根拠にはしない。
- J. C. Simo, T. J. R. Hughes (1998), *Computational Inelasticity*, Springer, DOI [10.1007/b98904](https://link.springer.com/book/10.1007/b98904)。出版社の書誌・目次（積分算法 pp.113–153、多面塑性 pp.198–218、一般 return mapping pp.219–239）のみ確認。本文は取得できなかった。return mapping と整合接線の本文を本書から引用しない。今回の具体的な診断根拠は上記 Sloan 本文と現行ソースに置く。
- R. I. Borja, S. R. Lee (1990), *Cam-Clay plasticity, Part I: Implicit integration of elasto-plastic constitutive relations*, CMAME 78,49–72, DOI [10.1016/0045-7825(90)90152-C](https://www.sciencedirect.com/science/article/pii/004578259090152C)、および R. I. Borja (1991), Part II, CMAME 88,225–240, DOI [10.1016/0045-7825(91)90256-6](https://www.sciencedirect.com/science/article/pii/0045782591902566)。出版社検索で書誌を確認したが本文取得は403。書誌のみ。今回: MCC の受理範囲をまず実CLIで確認し、未確認の算法を追加しない。
- P. Wriggers (2006), *Computational Contact Mechanics*, 2nd ed., Springer, DOI [10.1007/978-3-540-32609-0](https://link.springer.com/book/10.1007/978-3-540-32609-0)。出版社書誌と目次（接触構成則 pp.69–108、解法 pp.309–364）のみ。摩擦接触 return mapping の本文は取得できなかったため書誌のみ。T. A. Laursen, *Computational Contact and Impact Mechanics*, Springer, DOI [10.1007/978-3-662-04864-1](https://link.springer.com/book/10.1007/978-3-662-04864-1) も出版社書誌・目次のみ（同ページ copyright 2003）。両書から未確認の式や収束保証を引用しない。今回: 接触の離散更新と平衡経路の区別を現行コードで追う。
- Abaqus 2026 [UMAT](https://docs.software.vt.edu/abaqusv2026/English/SIMACAESUBRefMap/simasub-c-umat.htm) の Stability / Convergence Rate / PNEWDT と [Convergence Criteria for Nonlinear Problems](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEANLRefMap/simaanl-c-convergcriteria.htm) の材料計算失敗・平衡反復・severe discontinuity の節を本文確認。材料 Jacobian は全体 Newton の収束率に関わり、材料側からの切戻し要求と接触状態変化を含む全体反復の判定がある。材料局所積分を細かくするだけでは接触・全体経路誤差を制御したことにならない、という診断設計に使う。本コードに Abaqus の制御を移植したという意味ではない。
- Code_Aster U2.04.01 [§3](https://codeaster.gitlab.io/doc/docaster/manuals/man_u/u2/u2.04.01/Conseils_communs_aux_3_types_de_non-lin_arit_s.html) 本文確認。接触・材料・分解等の失敗による全体時刻列の再分割、Newton 残差、時刻列を2倍細かくして離散誤差を調べる方法が別項に記載。今回: 全体再求解の水準差と局所再生を区別し、失敗判定と残差を保存する。許容値緩和の根拠には用いない。

Code_Aster 同文書 [§6](https://codeaster.gitlab.io/doc/docaster/manuals/man_u/u2/u2.04.01/Non-lin_arit__de_contact-frottant.html) も本文確認。連成で問題が出る場合、各非線形を別々に確認するよう記載される。今回の片側固定・局所再生対照の動機とするが、追加拘束・接線変更等は行わない。

## 固定問題と観測行列

前回と同じ12 Tet4 / 16節点 / 12材料GP / 6接触GP、固定水位.75 m・γw=10、CONTINUOUS、非関連MC c=30 kPa、E=100000 kPa、ν=.3、φ=30°、ψ=5°、tension cutoffなし。open=.12、全7工程、各32増分を基準にする。物理入力・全体/線形残差の基準は変更しない。

材料上限m0/m1/m2=`1e-4/5e-5/2.5e-5`、接触上限c0/c1/c2/c3=`1e-4/5e-5/2.5e-5/1.25e-5 m`。m0〜m2×c0〜c2の9組合せ、m0c3/m2c3、m2c2の全体64、分割1対照を本体・独立で新規求解した。さらに、独立系で材料の実分割数8/32×接触の実分割数4/16を全増分・全GPで固定した4経路を解いた。**本体13＋独立17、計30経路・7,168受理状態**。全て完了、固定全体列からの切戻し0。ここで30は比較行列の経路数。別に観測なしwhole 1、元MC open=1.0のCLI/観測2、MCC3条件のCLI/観測/独立9を実行し、診断実験の新規完走求解呼出しは計42、受理増分は計10,040。2入力拒否・7実再開・局所再生・全件ソフトウェア試験中の求解はこの42へ加えない（`solve-call-ledger.json`）。MC最大力残差は本体9.90894061357e-9、独立9.61137431851e-11、真の線形残差は本体9.98337971992e-11で、元の1e-8/1e-10/1e-10を維持した。

独立は既存の入力組立3D MC系を固有名 `c03diag_*` に移した診断snapshot。幾何・水圧荷重・材料・接触から自由DOFの平衡を独自に解き、本体の状態・受理列・接線・実分割数を解の入力にしない。固定分割対照も各端点で全体平衡を求解する。これらは元の方式を置換する実装変更ではない。

旧封印の対応する6条件×2実装、計12経路と、今回の全履歴・全場JSONオブジェクトが**完全一致**した（`old-baseline-binding.json`）。旧保存値の再書換えや、新しい許容差による旧結果の置換はない。基準の変位/応力/塑性ひずみL2は元の3.438294562e-5 / 3.541744414e-5 / 1.505641374e-4を保持する。

## (i) 測り方とGP対応

元の `||a-b||₂/max(||b||₂,1e-12)`、全受理点・全成分の定義を保ち、絶対L2、分母、最大絶対差とその時刻/GP/成分、符号和、7端点、工程別値を併記した。反力・水荷重・応力・ひずみ・MC内部変数・接触の23場を保存する。ノルムは時刻→GP/節点→成分のC-order、Python `math.fsum` に固定し、Windowsの集計順差を閾値で吸収しない。新しい定義による受入は作らない。

本体の全13経路の材料/接触GP記述子は完全一致。材料はelement_index/local_index、接触はjoint_id/face_index/面内の固定順序で対応付ける。独立の位置差0、材料重みと接触の法線・重み・形状関数は既存の1e-12内で照合し、節点/所有要素対応も完全一致を検査した。最近接写像で別GPへ移し替えていない。

|接触隣接水準|本体のすべり相対L2|独立のすべり相対L2|
|---|---:|---:|
|c0→c1、材料m0固定|3.8839635e-05|4.513970017e-05|
|c1→c2、材料m0固定|0.0001114211284|7.900028122e-05|
|c2→c3、材料m0固定|1.300160876e-05|3.787703664e-05|

問題のc0→c1、c1→c2の**絶対L2**も本体1.36331610432e-6→3.91093148809e-6 m、独立1.58450200116e-6→2.77303976396e-6 mで増える。正規化の分母だけの見かけではない。端点と全履歴を混ぜない。元の7端点u=.005 / traction=.01ゲートは全13対応経路で通るが、全履歴・全内部変数・仕事の収束の合格へ拡張しない。全場の端点/全履歴/工程別値は `measurements.json`。

## (ii) 材料局所積分と接触局所積分

材料上限を1/4のm2へ固定しても、接触c0→c1→c2のすべり差は本体3.911271182e-5→1.119647703e-4、独立4.461290609e-5→7.969730199e-5と再増加する。m0からm2へ変えるとc1−c2差ベクトルは約2.01%（両実装）変わり、cos≈.9998。材料積分の寄与は存在するが、この対照では再増加の大部分を説明しない。

局所MCは、本体と同じ旧状態・ひずみ増分・実分割数を使う独立再評価で、前回同様に応力・塑性履歴が丸め程度まで一致した（`same-state-local-audit.json`）。加えて、生産MCカーネルを実分割の各サブステップで観測した。m0c0/m0c1/m0c2/m2c2本体の10,752点更新・**49,866局所サブステップ**はelastic 47,470、shear_face 2,396、角/頂点0。増分末の応力・塑性ひずみ・相当塑性ひずみ・塑性仕事・間隙比とカテゴリは本体と完全一致した（`return-trace-summary.json`、追加ZIP）。従ってこの再増加を角/頂点分岐の不一致に帰する根拠はない。未観測の一般経路について否定はしない。

保存された各変位/ひずみ経路を共通材料32分割で局所再生しても、c1−c2の応力差ベクトルの変化は約.232%、塑性ひずみ差は約1.15%。その応力・塑性差の大部分は保存経路の側に残る。一方、材料分割そのものは絶対場へ影響する。入力のみの独立全体平衡で片側を固定した値は次の通り。

|固定数の対照（独立全体平衡）|u L2|バルク塑性ひずみL2|JOINTすべりL2|
|---|---:|---:|---:|
|fixed-m8c4--fixed-m8c16|1.47467699e-05|4.08073477e-06|5.25931347e-05|
|fixed-m32c4--fixed-m32c16|1.47484589e-05|4.08028282e-06|5.26247881e-05|
|fixed-m8c4--fixed-m32c4|5.79629445e-05|0.000273297473|0.000237990436|
|fixed-m8c16--fixed-m32c16|5.79634938e-05|0.000273300538|0.000237851871|

片側凍結とは、材料/接触の分割数または保存入力履歴を固定する対照である。材料を別物理へ置換したり、接触係数を変更したりしていない。材料局所積分の影響がゼロとは結論しない。

## (iii) 接触局所更新と全体平衡の経路

同じ保存jump履歴に接触の共通32分割を使った局所再生では、問題のc1−c2すべり差は次のようになる。

|実装|元の差の絶対L2 m|共通32分割にした経路間差 m|元に対する残存ノルム|差ベクトルのcos|
|---|---:|---:|---:|---:|
|native|3.91093148809e-06|2.40327413832e-06|0.614501723|0.884876373|
|independent|2.77303976396e-06|1.75287775935e-06|0.632114181|0.886729123|

共通16→32分割によるこの経路間差の変化は元差の約.060% / .064%まで小さく、この固定された保存経路上の診断値は局所分割だけではほぼ変わらない。しかし元差の約61% / 63%が経路側に残り、元の局所分割を共通化することによっても差ベクトルが約54% / 53%変わる。**接触局所積分と全体で得られた経路の両方に寄与がある。** これらは相関を持つベクトルの比較で、足して100%になる因果割合ではない。局所再生は力平衡を解いておらず、参照解や全体平衡の代用にしない。

本体の分割数はNewton候補に依存して増大し、拘束開始状態へ戻って反復する。独立系は収束候補から必要数を増やす。観測した実分割ヒストグラムとiterate/repartition/line-searchは封印に保存した。分割1対照の本体対独立すべりL2は7.006498353e-9で、分割方針と端点間更新の違いを切り離す対照となる。ポリシーを完全に入れ替えた全体反実仮想ではないので、唯一の原因や一般保証とはしない。

同じm2c2の全体32→64を入力から解き直すと、共通時刻のすべりL2は本体.01011297754、独立.01013803643。接触の隣接内分割差より大きい。これも全体経路を細かくするだけで元の基準に合格したという意味ではなく、最細局所分割にも全体経路の離散化が残ることを示す。

## 元open=1.0の棄却箇所

旧2026-11-20等の切戻し/停止と失敗記録は保持する。今回はc=30の基準からopen係数だけを元の1.0へ戻した同じ小型入力を、実CLIと観測exampleで別々に解いた。**全工程まで完了したが、自動切戻し2回**。地盤初期化32点＋施工284点＝計316受理点で、固定32増分列の確認条件は未達。

|棄却工程|progress区間|停止判定|残差|局所GP情報|
|---|---|---|---:|---|
|open|.03125→.0625|Newton上限50|.00503398211736|材料GP6/7/8/11がshear_face。他はelastic。最大ひずみ増分GP6|
|open|.046875→.0625|全体line search失敗、iteration26|.00499626846951|同じ塑性GP。最大ひずみ増分GP8|

全体の判定は自由節点の縮約力残差で、**単一の「失敗GP」はない**（nullで保存）。最大縮約残差DOFは19/23。接触GP全6点のjump/gap・実分割数、全材料GPのreturnと分割数を `open1-stops.json` に保存した。最初の棄却では材料最大16,506分割、接触最大3,414分割へ拡大していた。材料局所の範囲拒否・評価errorは観測0。これは数値的な全体不収束による棄却であり、物理的破壊を証明したものではない。係数・正則化・許容値・拘束を変えて原因を消していない。

実CLIと観測exampleは実行版が異なるため、Run全体の完全一致検査ではbinary SHAの1項だけが異なった。そのfalseも保持し、全数値payloadは時計と実行版識別子以外完全一致と別に検査した。これを同一実行版での実再開完全一致と混同しない（`scope-payload-equivalence.json`）。

## MCCの現行受理範囲

既存MCCのBackward Euler、M=1.2、λ=.2、κ=.04、ν=.3、基準圧1 kPa、基準NCL比体積3.5、初期等方有効応力100 kPa・e=.7・塑性履歴0、pc=200/105 kPaを明示したbody_states入力を使った。水位・JOINT則・拘束・全体刻みは同じ。標準 `soil-fem-stages/1` の3ケースは実CLIで受理され、全7工程・224受理点・切戻し0。受理範囲だけ入力のみ独立3Dを新規求解した。

独立MCCは同じBackward Euler離散式のqとpcを消去して2未知数をSciPy HYBRで解く。生産の解析4×4 Newton・接線・解を呼ばない。新しい材料物理・接続は追加しない。現行MCC状態にない相当塑性ひずみ/塑性仕事を作らず、pcを含む実在する21場を全履歴比較した。

|MCCケース|塑性ひずみL2|pc L2|すべりL2|本体の状態|
|---|---:|---:|---:|---|
|mcc-pc200|0|0|0|open=.12、弾性、全1344 GP状態stick|
|mcc-pc105|0.01357368723|2.113370999e-06|0|open=.12、塑性あり、全1344 GP状態stick|
|mcc-pc105-open1|0.001434776687|2.490044688e-05|7.357607239e-05|open=1.0、塑性あり、slip106・stick1238 GP状態|

追加open=1.0ケースの本体最大塑性ひずみは.00365172718。この3入力では開口GP状態は0であり、MCCの開口履歴や一般化を確認したとはしない。MCCの別プロセス実checkpoint再開も今回の7地点MC再開へ含めない。CLI受理と全場収束を区別し、上の差を新しい許容値で合格へ変えない。

`soil-fem-mixed-mcc-stages/1` に同じ12 Tet4・水/JOINT入力を渡す組合せは実CLI/観測の両方で拒否された。理由は「mixed MCC requires dry Tet10/order>=3 ... water ... joints ... need explicit mixed integration」。拒否を保持し、新しい混合MCC接続は実装しない。一般のMCC水理/混合要素の受理を主張しない。

## (iv) 最細参照の評価

c3を厳密解にしない。材料m0の接触差は本体3.88396e-5→1.11421e-4→1.30016e-5、独立4.51397e-5→7.90003e-5→3.78770e-5。最後の差は縮むが、本体の隣接差cosは.216→.389、独立は.453→.978で、両実装が同じ漸近誤差場に達したとは確認できない。m2でも再増加が残る。

最細m2c3の本体対独立すべりL2は3.740388922e-5、バルク塑性ひずみは1.184517469e-4。一方、全体32→64のすべり差は約1%。最細局所条件の誤差上界、連続解、空間極限は未認定。Richardson補正・外挿・参照値の置換はしない。Tet10・長系列・96/768/6144のPCコマンド、予算、出力と見る値は [REPRODUCE.ja.md](REPRODUCE.ja.md)。PCは未実行。

## 実checkpoint再開、封印、判定

基準MCの各7工程progress=.5の実ファイルから別PIDで再開し、最後まで208/176/144/112/80/48/16、計784後続状態が完全一致。全Runもelapsed_seconds以外完全一致。同じ実行版SHAを要求し、観測なしwholeとの全Run一致も確認した。これは封印再監査と別の実行である。

新規求解・比較経路の全確定状態・Newton/分割/棄却・独立場・実行版・入力・checkpoint・再開・局所対照・使用ソースは `evidence.zip.part*` と `archive-index.json`、材料全サブステップの追加観測は `return-trace-evidence.zip` に保存した。全part/ZIP/全エントリーのSHAを検証し、封印ソースを隔離プロセスで再計測して `measurements.json` のJSONオブジェクト完全一致を確認した。封印確認を新規求解や実再開へ数えない。求解時の新exampleは後のrustfmt前、新studyは既定実行版名の編集前の原本も保存し、旧数値ソースは不変。最終試験のソース束は別に固定する。

|寄与|今回の根拠付き判断|
|---|---|
|測り方|絶対L2でも再増加。分母・端点だけ・GP誤対応だけでは説明しない。元定義を保持。|
|材料局所積分|絶対場には影響するが、m0→m2で再増加差ベクトルの変化は約2%。対象の全サブステップに角/頂点なし。同一状態・実分割のカーネル不一致は示されない。|
|接触局所積分|共通局所分割に変えると差ベクトルが約53〜54%変わる。寄与を限定確認。|
|全体経路|共通接触32分割でも元差の約61〜63%が残る。全体32→64で約1%のすべり差。候補依存分割と平衡履歴の影響を無視できない。|
|独立参照|独立側にも再増加があり、最細の全体刻み依存が残る。最細を真値と認定しない。|
|単一の主要原因|**判定不能**。接触局所積分と全体経路が主要な観測寄与で、材料積分・参照分割との連成を一つの因果割合へ分離できない。この同定失敗を failure-record.json に保持。|

**依頼された原因切り分けの必須クラウド診断は実施済み。優先するクラウド残件はなく、次はPCで元条件のTet10・長系列・大規模と独立参照を検証する。** 任意の小型追加水準は可能だが、それだけで未確認の一般化・MCC開口・最細参照を埋めたとはしない。元基準を全て満たしておらず、C-03を完了扱いにしない。

新フォルダのacceptance/measurements/validationは、固定コミット・SHA-256・Git blob・boolean/null pointerでP6-01補助台帳へ登録する。旧台帳原本・旧登録を再生成しない。新しい試験は通常core群に属し、重い再求解をCIで起動しない。文書はread_textの改行正規化を使い、実CRLFコピーを試験する。NCコード・設定、旧失敗・封印・参照値は変更しない。

全件の終了コード・成功/失敗/未実行/skipと生ログは `validation.json` / `validation/`。最終検証値は同ファイルの確定結果を参照。push後のWindows CIは未確定。

## Linux全件検証の確定結果

fmt・全target clippy（警告をエラー）・locked release buildは成功。**Rust466成功 / 0失敗 / 0未実行 / 0skip、Python769成功 / 0failure / 0error / 0未実行 / 0skip**。Rustは標準463＋exampleの3試験、docは0件。Rustの発見/実行466 IDとPython全7群の発見/実行769 IDは重複・欠落なく一致した。expected failure / unexpected success / ignored / filtered-outも0。Pythonは通常340/405、専用4/4/4/4/8件で、新規8試験を含む。1270ソースは全Python実行の前後でSHA一致し、163 Rustソースは試験snapshotと最終checkoutで全件一致した。

初期の書式・PowerShell書込先・大容量Git取得のHTTP/2切断・Gitオブジェクト復旧/履歴取得の準備上の不成功は、最終成功件数へ加えず `failure-record.json` と `validation/preparation/` に保持した。skipで回避せず必要な証拠を復旧し、厳格なP6事前監査は301 source・integrity_ok=true。最終検証snapshotを固定した後の台帳監査と関連再確認は別ログに置き、769件へ二重加算しない。push後の新先端Windows CIは未確定。

台帳登録後の関連再確認は91件すべて成功（失敗・未実行・skip 0）。P6最終監査は301 source・integrity_ok=true。N047/N048は科学証拠コミット `2f25ea3f68c9375ae67941b15933ae2b415ccf94`、N049は全件検証コミット `059ef2e16453b34bc06391e5f048c4b69251005c` を固定し、登録後の結果は `validation/post-registration.json` と `validation/p601-final-audit.json` に別記した。旧7756ファイル中、変更は台帳の追加登録と監査SHA定数の2ファイルだけで、残る7754のGit blobは完全一致する。

C-03は未達、passed=false。
