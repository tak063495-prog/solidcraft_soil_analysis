# A-01：乾燥・複数JOINT／軸力部材複合の独立Python参照

起点は PR #13 の `codex/plastic-spatial-convergence`、指定コミット `7ae0a2a669e0e1df93f3dd9791bbeda056334c3e`。専用ブランチ `codex/py-multi-joint-structure-composite-reference` への新規追加だけである。保存先の日付ラベルは依頼指定、実施日は2026-11-12（日本時間）。**A-01・NC-01・P6は未達を維持し、拡張範囲の `accepted=null` を保持する。draft PRをマージしない。**

## 範囲と既存記録の確認

[元の受入表](../2026-11-13-joint-ext-combos/README.ja.md)、[acceptance.json](../2026-11-13-joint-ext-combos/acceptance.json)、[A-01](../2026-11-16-issue-triage/todo.ja.md) を読む。`scope_audit.py` は元の `evidence.zip` 自体と範囲拒否48エントリーを既存manifestのSHA-256で照合した。保存された12拒否ログと調査対象ソースのSHAは [scope-audit.json](scope-audit.json)。本体・封印された実行ファイルは、この調査・新参照から実行していない。

|対象|確認できた条件|根拠・限界|
|---|---|---|
|既存の限定受入|乾燥、単一の先端宣言付き適合Tri6面集合、全所有要素活動、A=`transfer_force`／B=`stress_free`、接合専用DEATH→再誕生|元受入表とcombination_1。`passed=true` はこの限定範囲だけ。`requested_all_combinations_passed=false` を保持|
|複数JOINTの現行拒否条件|いずれかに `birth_policy`、いずれかに `crack_tip` があり、`joints.len()!=1` または先頭JOINTにwaterがある|`src/joints.rs` の “declared-tip JOINT BIRTH is limited to one dry interface”。ソース読取で確認した条件。新しいCLI拒否実行とはしない|
|構造複合の現行拒否条件|同一JOINTに `crack_tip` と `birth_policy` があり、`beam_system` または `anchor_system` が存在|`src/stages.rs` の “declared-tip JOINT BIRTH with structural members is unverified”。`src/beams.rs::Configuration` は plates も収容する。板だけの具体的な旧入力・ログは未確認|
|保存拒否入力|`evidence.zip/scope-rejections/` は湿潤BIRTH9件とRichards3件|全12 caseを点検し、複数JOINTまたは構造systemを含む対象入力は0件。**乾燥の複数JOINT／梁・板・アンカー複合の正確な旧入力・CLI拒否ログは未確認**。範囲外宣言やソースガードを実行ログに読み替えない|
|一般の乾燥JOINT|先端宣言／BIRTHを持たない複数JOINTまで一律拒否するとは確認していない|上記ガードの適用範囲を一般化しない|

旧 `retained-failures/`、限定受入、参照値、許容値を保存する。新参照の固定世代接触は、未受入の先端BIRTHや初期力投影を受理する実装ではない。

|並行作業|今回の境界|
|---|---|
|#16 受入試験|既存収集・CIの件数条件を変更せず、strict xfail/skipを追加しない|
|#17 湿潤JOINT、#19 過渡水理、#21 湿潤土BIRTH、PR #13 C-04|waterを持たない乾燥だけ。水圧・水位・Darcy・Richards・水BIRTHを実装しない|
|#18 Griffiths、#24 MCC|均質等方弾性bulkだけ。斜面、SRM、MCC・MC塑性bulkの参照へ拡張しない|
|#20 曲面殻、#22 平面Tri6オフセット複合、#25 折れ目/段差|板・殻・オフセット・director・折れ目幾何を求解しない。板複合は未確認のまま|
|#23 単一JOINT先端BIRTH投影|世代0固定、BIRTH/DEATH・初期力投影・先端力の回収を実装しない。共有節点は通常の1物理節点|

## 独立実装と出典

`reference.py` は標準ライブラリ、NumPy、SciPyだけを使う。Rust、既存Pythonソルバ、既存Runの数値を呼び出さず、自作Tri3平面ひずみ／Tet4、小ひずみ弾性、乾燥接触、線形軸力部材を組み立てる。`audit.py` は求解・assemble・接触return関数を呼ばず、入力幾何と生変位からテンソル応力の発散、接触式、軸力、仕事を照合する。既存Python全件試験が本体を実行することと、新参照の独立性を区別する。

|式・規約|出典と照合|
|---|---|
|`σ=λ tr(ε)I+2Gε`, `λ=Eν/((1+ν)(1−2ν))`, `G=E/(2(1+ν))`|[GetFEM公式・等方線形弾性](https://getfem.org/userdoc/model_linear_elasticity.html)。小ひずみ・平面ひずみの定義を確認。自作のaffine patchと剛体並進で検査|
|Coulomb stick/slip、penalty、trial/return、後退Eulerの履歴更新|[COMSOL公式・Tangential Contact with Friction](https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.088.html)。一般原理の出典であり、COMSOL実装全体との一致を主張しない。下記の式を本参照のモデル定義として自前実装|
|Gauss–Legendre節点と重み|[NumPy公式 leggauss](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html)。線分を[0,1]へ写像。Tri3は下記Duffy変数変換を自前導出し、面積・一次モーメントを別検査|
|軸力と仕事|`δ=(u_b−u_a)·a`, `N=EAδ/L`, 端内力 `(-Na,+Na)`, `U=Nδ/2`。一次元の `ε=δ/L`, `σ=Eε` と体積ALから直接導出。自作の75 kN閉形式試験で照合。梁の曲げ・ねじり式は使わない|
|既存診断基準|[元README](../2026-11-13-joint-ext-combos/README.ja.md) と acceptance.json の `tolerances`。本拡張の受入基準へ転用しない|

接触の法線nはminus→plus、jumpは `j=u_plus−u_minus`、`g=n·j` は開離正。出力traction qはplus側の**内力**で、実際にplus側へ作用する接触外力は−q。`P=I−n⊗n`, `j_t=Pj`, `p_n=−k_n min(g,0)`。閉鎖時のtrialは `r=k_t(j_t−s_old)`。`||r||<=μp_n` なら固着、超える場合は `q_t=μp_n r/||r||`, `s_new=j_t−q_t/k_t` とする。`q=−p_n n+q_t`。開離はq=0、旧塑性すべりを保持する。これは明示したCoulomb構成則のreturnであり、生結果の平滑化・クリップ・補正ではない。Newton反復中に履歴を確定せず、平衡後だけcommitする。

bulkの節点内力は `∫Bᵀσ dV`、接触内力は `±∑wNq`。別監査では `V σ grad(N)` を直接散布する。ペナルティに伴う弾性蓄積も `U_c=∑w(k_n min(g,0)²+||q_t||²/k_t)/2` として含める。有限ペナルティ剛性のモデルであり、非貫入拘束の連続極限の検証ではない。

Tri3の積分は `N=((1−a)(1−b), a, (1−a)b)`, `0<=a,b<=1`, 面積Jacobian `||(x1−x0)×(x2−x0)||(1−a)`。1次のDuffy積1点は重心則と同一ではない。GP系列は積分次数の変更で、メッシュh系列や点ごとの空間収束に読み替えない。

仕事は生の反力を用いた台形積分 `W_ext=Σ(f_old+f_new)·Δu/2` と、`D_f=Σw μp_n,new ||Δs||` をそれぞれ積算する。**`W_ext−(U_bulk+U_c+U_member)−D_f` は補正せず保存**する。閉鎖両端点では、接線ばね関係から台形仕事の離散恒等式

`W_contact,trap − ΔU_c − D_f + (w/2)(q_t,new−q_t,old)·Δs = 0`

を直接導出できる。これを独立照合するが、最後の項を外力仕事へ足して物理的仕事合格にすることはない。開離端点にこの恒等式は適用せず件数を保存する。時間離散の経路差、ペナルティモデル誤差、GP積分差を分ける。

## 10ケースと負荷

E=3000 kPa、ν=.25、k_n=50000 kPa/m、k_t=30000 kPa/m、μ=.25、2D厚さ1 m。1 m角のブロックを領域内で接続し、接触両側の節点を複製する。parallelは3列（JOINT2本）、cross/sharedは2×2、Tは下部2 m幅＋上部2ブロック。3Dは同じ配置を奥行き1 mへ押出す。sharedだけ中心点（3Dは中心縦辺）を1物理節点にする。接触面が交差するだけのcrossは中心の各ボディ節点を分離する。

外周節点の全並進を規定し、内部・界面節点は自由。初めに `u=−.004x` へ圧縮、次に `u_x+=.032 y f`, `u_y+=.022 x f`（3Dは `u_z+=.013 x f`）でf:0→1、最後にf:1→−.5へ反転する。各工程8/16/32増分。位相ごとに連続な経路で、共通の実端点を比較し、補間で一致させない。

beam2d/anchor3dは接触節点に2本の線形軸力部材を結ぶ。EA=600/300 kN。アンカーも引張専用化せず、両方向の線形軸力モデルである。1本は界面内部、もう1本は近傍の外周節点にあり、近傍への伝達はbulk要素そのものが担う。最近傍snap、隠れた剛体リンク、梁曲げ・板オフセットは使わない。

## 小型測定

80 Run、6,480の確定増分を保存した。基準は各工程32増分・GP次数2、bulk分割2/m。以下は基準Runの全履歴最大値。全Runの生値と入力メタデータは `evidence.zip`、SHAは `results/study.json`。各Runの `complete=true` は独立参照計算の完了でありA-01受入ではない。

|参照|節点 / 接触GP|最大節点残差 kN|最大外力モーメント残差 kN m|終点の仕事差 kN m|
|---|---:|---:|---:|---:|
|parallel2d|27 / 8|5.773160e-14|1.634248e-13|-7.879441e-04|
|cross2d|36 / 16|1.358913e-13|2.753353e-13|-2.964498e-03|
|t2d|33 / 12|1.376677e-13|1.740830e-13|-1.513945e-03|
|shared2d|33 / 16|1.243450e-13|2.309264e-13|-1.179625e-03|
|parallel3d|81 / 64|1.210421e-11|4.298784e-13|-3.100359e-04|
|cross3d|108 / 128|1.928924e-11|2.913225e-13|-4.824052e-04|
|beam2d|27 / 8|5.240253e-14|2.131628e-13|-8.767343e-04|
|anchor3d|81 / 64|7.575107e-12|2.557954e-13|-2.376794e-04|
|t3d|99 / 96|1.441591e-11|2.078338e-13|-2.399626e-04|
|shared3d|99 / 128|1.829106e-11|2.282619e-13|-4.566917e-04|

別監査の最大差は、テンソル応力からの全節点内力1.421086e−14 kN、接触力散布0 kN、軸力0 kN、Coulomb式3.552714e−15 kPa、接線ばね関係1.735570e−18 m。摩擦仕事恒等式3.509188e−19 kN m、閉鎖両端点の離散仕事恒等式1.216736e−17 kN m。前者の再組立と後者の代数恒等式は既存1e−10の診断内である。**物理的な外力仕事の差は残り、全すべり履歴を既存の絶対1e−10仕事基準で合格にしない。** モーメント・すべり・順序の拡張受入許容値は未定義。

beam2dの右ボディ終点x成分は、接触内力−14.8046395045953 kN＋部材内力−4.66705460264413 kN＝反力−19.4716941072394 kN（bulk合力は丸め差）。anchor3dでは−10.5258530083659−2.48336897397982＝−13.0092219823457 kN。y/z成分・部材端力・ゼロ端モーメント・各エネルギーも全節点台帳に残す。この限定した軸力分担を一般梁／板の合格へ転用しない。

### 評価順・JOINT順・節点番号

接触GPの逆順、JOINTブロックの逆順、seed=20261121の節点番号置換を**それぞれ新規求解**し、物理IDへ対応させて全履歴を比較した。10ケース×3変更の世代・open/stick/slip状態はすべて一致。全変更の最大差は変位8.734333e−15 m、塑性すべり8.456343e−15 m、traction 1.231353e−10 kPa、反力3.383072e−12 kN、外力仕事3.996803e−15 kN m。数値はビット完全一致ではなく、生差を `results/observations.json` に保存した。新しい順序許容値は未定義、`accepted=null`。単なるCSV行ソートの試験と、この再求解を区別する。

### 荷重刻み・GPの観測系列

|参照|終点仕事差 8→16→32増分（kN m）|GP次数1→2 / 2→3の全履歴節点変位相対L2差|
|---|---|---|
|parallel2d|-6.198143e-03 → -2.263177e-03 → -7.879441e-04|1.231232e-03 / 6.124768e-17|
|cross2d|-3.058148e-02 → -7.777571e-03 → -2.964498e-03|1.940070e-03 / 1.164616e-16|
|t2d|-7.202388e-03 → -2.722751e-03 → -1.513945e-03|8.260364e-03 / 6.706287e-04|
|shared2d|-1.113906e-02 → -2.888663e-03 → -1.179625e-03|1.598010e-03 / 9.455912e-17|
|parallel3d|-2.180402e-03 → -7.682717e-04 → -3.100359e-04|9.203631e-04 / 4.079038e-15|
|cross3d|-5.035358e-03 → -1.962952e-03 → -4.824052e-04|8.491980e-04 / 2.992206e-06|
|beam2d|-5.614395e-03 → -2.176728e-03 → -8.767343e-04|9.980659e-04 / 6.149649e-17|
|anchor3d|-2.462048e-03 → -6.003873e-04 → -2.376794e-04|7.800757e-04 / 2.573892e-15|
|t3d|-3.706740e-03 → -1.003112e-03 → -2.399626e-04|1.895122e-03 / 8.353094e-05|
|shared3d|-3.329289e-03 → -1.349550e-03 → -4.566917e-04|1.268194e-03 / 3.234752e-15|

仕事差の大きさはこの3水準で減少するが0に到達せず、合格値や連続解を作らない。全履歴のすべり相対L2差は、T2Dで0.0399190→0.0342439（記述的次数0.221229）、T3Dで0.0478922→0.0170993（1.485852）。3D平行は0.000527201→0.000402360、3D交差は0.000308717→0.000193975。全フィールド・負次数もJSONに保存し、1e−16程度の丸め差から得た次数は収束の証拠にしない。

GPの非一致点は補間しない。GP系列では全節点場・部材軸力・JOINTごとの明示した面積積分・仕事だけを比べる。T2D/T3D/交差3DのGP2→3差は残る。これは次数による観測差であり、点ごとの接触場収束、h収束、先端特異場、実CADでの合格は未確認。Richardson補正・外挿参照は生成しない。

## 将来の本体CSV/JSON比較

`compare.py` は本体実行を行わず、正規化された生出力だけを読む。契約の例は `evidence.zip` 内の完全なJSON。schema、物理・離散化、単位、符号、全節点/部材/JOINT/GP ID、接触世代、開閉・すべり、接触力、断面力、全節点力と仕事の全履歴を必須にする。JOINT IDの対応は明示の全単射JSONだけを認める。欠損・重複・未対応ID、非有限、途中状態、無申告のraw、平滑化/クリップ等のtransformations、単位/符号の自動推測、補間、GPの最近傍対応を拒否する。

CSVは `step,kind,identity,payload_json` の4列と `.manifest.json`。vectorと単位を列名から推測せず、JSON payloadの全成分を保存する。manifestにはフレームのスカラー台帳、全ID、評価順、CSVのSHA-256を持つ。JSON→CSV→JSONの完全往復を試験する。

比較器は**順序を変えた別出力を必須**とし、全生差を返す。新しい接触すべり・開閉・断面力・順序の比較許容値は未定義のため、数値差が1つでもあれば `RAW_DIFFERENCE_REQUIRES_DEFINED_POLICY`／exit 2として比較受入を拒否する。丸め差を含む実際の順序変更出力もこの厳しい契約では受理しない。全生値が一致する人工往復だけはデータ契約検証exit 0、数値受入は常にnull。既存の再組立1e−10 kN、仕事1e−10 kN m、独立弾性相対1e−7は出典付き診断に限り、未定義の比較を受理するfallbackに使わない。

最終値だけから、虚偽のraw宣言を伴う**すべての非公開補正を識別することは不可能**である。差・力台帳不整合・補正宣言・provenance欠落は拒否できるが、本物の本体生export・入力指紋・独立再組立が必要。現時点の本体exporter/adapter、実CAD入力との一致は未確認。本体と数値比較したとは主張しない。

## 試験・保存・未達

新規試験は `checks_reference.py` の明示実行だけで、既存の `tests/`・`tools/`、pytest設定、CI・件数条件を変更しない。skip/xfailで通さない。最終の件数・終了コード・旧証拠不足等の初回エラーは `results/validation.json` と `logs/` を参照する。初回の成功件数を再実行へ加算しない。新しい大規模3D・実CADはクラウドで実行せず、[Windows手順](REPRODUCE.ja.md)へ分けた。Rust全試験は今回の範囲外、既存Pythonのためのlocked release buildと既存Python全件を区別して記録する。

A-01の先端BIRTH・水圧・世代変更・独立一般3D全場・本体との比較・実再開は受入未達。板複合、一般梁曲げ、実CAD、大規模h/GP/刻み系列も未確認。NC-01、P6、#17/#21/C-04等の未達・過去失敗を完了に変更しない。新先端Windows CIは未確認で、Linux試験や旧CI成功、時間切れ・キャンセルを成功へ転用しない。

### Linux最終実行結果

[validation.json](results/validation.json) に全ID・割当・実行・終了コードを保存した。新規は**26件成功**（subtestを加算しない）、失敗・エラー・skip・xfail 0。既存は**543件を発見し4群すべての実行を試み、通常群272＋263＝535件が成功**。残る2証拠群はsetUpClassで各1エラーとなり、計8メソッドは未実行。**全件成功とはしない**。skip/expected-failure/unexpected-successは0で、試験やCI件数条件を変更して通していない。初回の通常群263/253実行は不足旧証拠とCargo PATHのエラーとして保持し、最終件数へ加算しない。

- 平均応力証拠群：封印archive・生JSON再構築の検査後、既存の完全一致条件 `saved mean-stress/complement partition differs from own history` で停止。SHAで照合した旧ソース・controlsだけの追加診断では29分配に差があり、最大7.105427357601002e−15 kN。通常/Haswell/Sandybridge/NehalemのBLAS指定でこの差は解消しなかった。新しいNC求解・力配分調整は行っていない。
- 空間証拠群：archive SHAと離散構造の検査後、元JSONとのバイト一致で停止。独立に保持した再計測JSONの差は浮動小数点49値、最大4.440892098500626e−16。非浮動小数点の差は0。元値や許容値を書き換えず、Linuxで `allow_platform_float_differences` を有効にしない。

環境起因の可能性はあるが原因は**未確認**。失敗・旧値・生差は `retained-failures/` に保存し、これら8件を受入や成功に数えない。locked release buildは成功、Rust全試験と新先端Windows CIは未確認。既存5,707パスのmode/blobを保持し、配置済み4,893ファイル（7,951,222,450 bytes）を元Git blobに照合して差0。疎チェックアウトの未配置814ファイルは再検証済みとはしない。既存ソースの終了時SHA一致は [final-preservation.json](results/final-preservation.json)。
