# 杭だけの独立 Python 参照と診断（P6-03 の未達を保持）

指定起点は `7ae0a2a669e0e1df93f3dd9791bbeda056334c3e`（PR #13、`codex/plastic-spatial-convergence`）。保存先の日付は依頼された整理ラベルであり、その日の実行を主張しない。既存ソース・テスト・CI・P6進捗表・既存 interim フォルダーを変更せず、このフォルダーの新規ファイルだけを追加する。

**杭の空間収束・先端/周面分担・塑性/すべり・杭列/群杭・P6は未達のまま。** A-01、NC-01、C-04にも受入を追加しない。MCC（#24）、水理/降雨（#19）、湿潤JOINTと他の並行参照PR #16〜#25 のモデルには触れない。新しい参照の工学的許容値は既存文書から特定できず **「未定義」、`accepted=null`**。既存の限定受入と今回の検算法確認を、適用範囲拡大の合格へ転用しない。

## 元のケースと判定の読み取り

[ToDo P6-03/P2](../2026-11-16-issue-triage/todo.ja.md)、[P6_STATUS P3b-05/P6-03](../../docs/P6_STATUS.ja.md)、[残件表 P3b-05](../../docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md)、[SOLID_PILE「固定した物理条件」「残る検証」](../../docs/SOLID_PILE.ja.md)を照合した。P6-03/P2の複数残件から杭だけを切り出す。P6-03のMCC性能記録を杭の性能実績と読まない。

|条件|確認内容・一次入力|
|---|---|
|形状・座標|地盤 `[0,10]×[0,10]×[0,8]` m。杭 `[4.5,5.5]×[4.5,5.5]×[2,8]` m、中心(5,5)、1 m角、長さ6 m。地盤794 m³、杭6 m³。[solid-pile.json](../../examples/solid-pile.json) のCAD長さはmm|
|材料・単位|m-kN-kPa。地盤 E=30,000 kPa、ν=0.3、γ=20 kN/m³、杭 E=30,000,000 kPa、ν=0.2、γ=24 kN/m³。[ground-solid-pile.json](../../examples/ground-solid-pile.json)|
|構成則の限定|元入力MCの弾性枝。c=1,000,000 kPa、φ=ψ=0、引張カットオフなし。旧限定検証は全GP塑性ひずみ/仕事0。今回の自作線形弾性はこの範囲だけに対応し、塑性/支持力を計算しない|
|界面・境界|共有節点による完全付着、周面と先端に独立のすべり則なし。側面法線ローラー（x側面 ux=0、y側面 uy=0）、底 z=0で uz=0。その他の表面は無載荷時自由|
|水・初期状態|乾燥、水位なし、初期応力0。杭が既に存在する全領域の自重平衡後を施工基準。K0施工履歴や打設履歴ではない|
|外力・比較工程|杭頭z=8、面積1 m²に(10,0,−100) kPa、合力(+10,0,−100) kN。自重平衡→載荷→除荷。[stages-solid-pile.json](../../examples/stages-solid-pile.json) と[ボディ平面版](../../examples/stages-solid-pile-body-planes.json)|
|比較量|杭頭面積平均の載荷増分沈下 `−1000*mean_A(u_loaded,z−u_gravity,z)` mm。節点単純平均や自重込み総沈下と混ぜない。水平も自重からの増分|
|旧離散化|436/819/2279 Tet10。公開wrapperのメッシュ指定は3/2、2/1.4 m。非一様Gmsh系列なので要素数から一様hや観測次数を推定しない|

[旧first-pair](../2026-10-26-body-pile/reports/solid-pile-first-pair.json) と[旧refined-pair](../2026-10-26-body-pile/reports/solid-pile-refined-pair.json) の数値は変更しない。

|Tet10数|元の追加沈下 mm|隣接する細側を分母とする差|
|---:|---:|---:|
|436|0.32521742815831706|436→819: 0.06687358294838157 = 6.687358294838157%|
|819|0.3485245109509388|819→2279: 0.04717393008743999 = 4.717393008743999%|
|2279|0.36577978075570766|さらに細かい旧結果はこの記録では未確認|

定義は `abs(s_fine−s_coarse)/abs(s_fine)`。**6.687→4.717%は連続体厳密解への誤差ではない。** 元文書は「空間収束は受入済みにしない」、JSONは `spatial_convergence_accepted=false`。新規の杭空間収束用数値しきい値は上記文書・公開入力/記録からは特定できず「未定義」。他ケースの5%などを持ち込まない。

旧限定受入は同一Tet10メッシュの独立3D全場・力/モーメント・VTU/生GP・小さいメッシュの再開など。[旧検証器](../../tools/verify_solid_pile.py) にある `relative<1e-7`（全節点変位）、`stress_error<1e-7`、`void_error<1e-8`、真残差 `<1e-10`、反力/全界面合力 `atol=1e-5 kN, rtol=1e-9` は**その同一メッシュ照合の定義**であり、空間収束しきい値ではない。今回のHex8/1D比較へ転用しない。今回この旧検証器は参照生成にimportも実行もしない。

## 独立参照・出典・未確認の対応

全スクリプトは標準ライブラリ、NumPy/SciPyと本フォルダー内モジュールだけを使う。本体Rust、`tools/fem_reference.py`、既存ソルバ/検証器、外部FEMパッケージは参照生成に使わない。出典は式を確認したものだけを次のように使用する。出典の実装コードの移植ではなく、以下の式・弱形式を自作で離散化した。

1. **軸力杭の閉形式**：Mosher & Dawkins, *Theoretical Manual for Pile Foundations*, ERDC/ITL TR-00-5 (2000), Appendix A, A1〜A6（A1〜A3頁、PDF物理頁105〜107）。[公開報告書PDF](https://vulcanhammer.info/wp-content/uploads/2017/08/erdc-itl-tr-00-5.pdf)。一定周面剛性の軸力杭と弾性先端条件を確認。記号を `EA u''−k_a u=0, N=−EA u', N(0)=P, N(L)=K_t u(L)` に統一し、`a=sqrt(k_a/EA), r=K_t/(EA a), t=a(L−s)` とおく。`u=P(cosh(t)+r sinh(t))/(EA a sinh(aL)+K_t cosh(aL))`、`N=−EA u'`、地盤から杭への周面反力/単位長さは `−k_a u`。これは出典の一定係数一般解に明示境界条件を代入した形。コードでは頭/先端条件、既知指数解、積分力釣合いを別途照合する。`k_a` はkN/m²、`K_t` はkN/mであり、摩擦応力kPaとの変換には周長が必要。
2. **水平Winkler梁**：Więckowski & Świątkiewicz, *Materials* 14 (2021), 460, §2.1, eqs.(1)〜(6), [DOI](https://doi.org/10.3390/ma14020460)、[確認したPDF](https://mdpi-res.com/d_attachment/materials/materials-14-00460/article_deploy/materials-14-00460.pdf)。`EI w''''+k_h w=0, M=−EI w'', V=−EI w''', r_h=−k_h w`。有限長、頭 M=0・V=−H、先端 M=V=0の定係数系を無次元座標 `t=s/L` の4×4行列指数で解く。形状関数を使わない半解析解と、三次Hermite梁＋分布ばねのGalerkin解を比較する。
3. **半無限解の範囲確認**：FHWA-HRT-06-106 (2009), ch.6 eqs.(4)〜(8), [FHWA原報告](https://www.fhwa.dot.gov/publications/research/infrastructure/pavements/06106/chap6.cfm)。式(8)の端モーメント0の場合、`beta=(k_h/(4EI))^(1/4)`、`w=H exp(−beta*s) cos(beta*s)/(2EI beta^3)`。導関数からM,Vを計算する。同頁式(7)のテキスト表示と式(8)には不整合があるため式(7)は使わない。式(8)を微分し自由端M=0,V=−Hと減衰を検算する。この報告はダウエルの梁式の出典として使い、その支持剛性を杭へ転用しない。半無限解を有限6 m杭の正解としない。
4. **同じ3D物理条件の小規模数値参照**：`solid_pile_reference.py`。線形等方弾性の構成則/弱形式は Logg & Langtangen, [The equations of linear elasticity](https://github.com/jorgensd/dolfinx-tutorial/blob/main/chapter2/linearelasticity.md) のPDE/variational formulation、Hex8形状関数は [FrontISTR文書「Shape Functions for First-Order Elements」](https://frontistr-documents-i18n-071b7b.gitlab.io/en/theory/discretization/08_element_solid.html) の一次六面体式を確認。`sigma=lambda tr(eps)I+2mu eps`、`K_e=∫B^T D B dV`、`N_i=Π(1+sign_i*xi)/8`。矩形要素の2×2×2 Gauss完全積分、材料界面で整合する共有節点、上記の自重/面荷重/ローラーを自作で組む。FEniCS/FrontISTRは呼ばない。弾性であるため自重解と載荷増分を別右辺で解き、和が載荷状態、増分を除くと自重状態に戻る。非線形載除荷履歴や保存再開の受入ではない。

1D参照の数値的適用範囲は `aL<=20`、有限梁の `(k_h*L^4/EI)^(1/4)<=5`。これを越える転送行列/双曲線関数の条件数対策は未実装で、明示エラーにする。大きなばね剛性を勝手に切り詰めない。

1Dの梁/棒要素の形状関数はそれぞれ節点値/勾配を満たす三次Hermite多項式と一次多項式を明示生成する。`∫EA u'v'+k_a uv`（先端 `K_t u(L)v(L)`）、`∫EI w''v''+k_h wv` を4点Gaussで完全積分。誤差の積分は要素ごとに独立6点Gauss。せん断力は平滑化せず三階導関数から計算する。数値行列の真残差も記録する。

**3D→1D対応の未確認値**：有限10×10×8 m地盤・正方形完全付着杭に対する `k_a, K_t, k_h`、外側/底面境界と地盤相互作用を含む換算の出典は未確認。杭のEA=30,000,000 kN、EI=2,500,000 kN·m²（正方形I=1/12 m⁴）、L=6 m、P=100 kN、H=10 kNまでは確認できるが、地盤Eをそのままばねkにしない。等面積円への置換や半無限地盤化も行わない。元3Dケースの閉形式値、梁杭との同条件一致、先端/周面分担は未達。この不足を自動同定や過去沈下へのfitで埋めない。

## 実行結果と位置づけ

[1D検算法結果](results/reference-results.json) は **明示的な検算法確認専用条件** `aL=1, K_t=EA/L, beta*L=1` を使う（元3Dの同定値ではない）。L/EA/EI/P/Hは上記、`k_a=EA/L²=833333.3333333334 kN/m²`, `K_t=5000000 kN/m`, `k_h=4EI/L⁴=7716.049382716049 kN/m²`。全行 `accepted=null`。

|n|h m|軸変位 相対L2|軸力 相対L2|梁変位 相対L2|梁M 相対L2|梁V 相対L2|
|---:|---:|---:|---:|---:|---:|---:|
|4|1.5|3.285961e-3|7.195871e-2|1.448192e-5|4.700731e-2|3.899479e-1|
|8|0.75|8.239592e-4|3.605807e-2|8.927942e-7|1.195019e-2|1.976885e-1|
|16|0.375|2.061448e-4|1.803890e-2|5.557680e-8|2.999813e-3|9.918059e-2|
|32|0.1875|5.154589e-5|9.020686e-3|3.538418e-9|7.507178e-4|4.963224e-2|

`p=log(e_h/e_h/2)/log(2)`。最後の観測次数：軸変位/周面反力 **1.999729**、軸力 **0.999803**。梁変位/地盤反力 **3.973306**、回転 **2.999647**、M **1.998530**、V **0.998780**。梁の変位誤差最細は丸め/条件数の影響を含むので理想4次への強制補正をしない。閉形式の軸杭頭沈下は0.020000 mm、有限梁頭変位は0.872172807 mm、半無限は0.432000 mm。有限/半無限の差は離散化誤差ではなく端条件のモデル差。

[同物理条件の小規模3D結果](results/small-3d-results.json)：x/yの分割基点0,4.5,5.5,10、zの分割基点0,2,8の各区間をm等分する入れ子Hex8系列。各m倍増で全方向h半減。h最大は最大辺長（対角長ではない）。

|m|Hex8数|節点|自由DOF|最大辺h m|杭断面h m|追加沈下 mm|
|---:|---:|---:|---:|---:|---:|---:|
|1|18|48|80|6|1|0.199839800511|
|2|144|245|546|3|0.5|0.290057891742|
|4|1152|1521|3926|1.5|0.25|0.346173412600|

`log2((s_1−s_2)/(s_2−s_4))=0.6850169053`。三点の沈下差による次数であり、厳密誤差に基づく次数ではない。漸近域・連続体参照の精度は未確認。外挿値/合格値を作らない。各水準の地盤/杭体積794/6 m³、頭面積1 m²、自重合力−16024 kN、載荷増分(+10,0,−100) kNを照合。最大自由DOF相対残差は約2.54e-12。実測約1.44秒の小規模計算であり、10⁴/10⁵級3D系列はクラウドで実施していない。

[旧記録の診断](results/historical-diagnosis.json) は旧ファイルのSHAと元比率を再読込して保持する。

|独立参照|旧436/819/2279の位置づけ|誤差源・制限|
|---|---|---|
|軸力杭閉形式/1D FEM|全て**判定不能**|元の周面/先端剛性と有限3Dへの対応が未確認。0.02 mmを元の沈下と比較して合否を付けない|
|有限/半無限Winkler梁|全て**判定不能**|水平の参照であり鉛直沈下と比較しない。水平ばね/有限境界への対応も未確認|
|同条件Hex8小系列|工学的判断は全て**判定不能**|同じ弾性連続体条件だが、近似離散解で誤差包絡/厳密参照ではない。Tet10との補間次数・非一様メッシュ・界面/角部勾配の違いが残る|
|Hex8計算済み値の単なるmin/max範囲|436値は**含まれる**、819/2279値は**外れる**|この範囲は誤差帯ではない。「外れる」を物理モデルの棄却や過去限定受入の失効と解釈しない|

公開した旧2件はsummaryで、全節点変位/全GP場は同梱されていない。旧同一メッシュ全場合格という記録を保持する一方、今回旧本体を再計算・再監査したとはしない。Hex8最細との差も厳密誤差ではない。先端/周面の共有稜線の節点力を恣意的に分けない。支持力・すべり・群杭/杭列換算・大規模性能は残る。

## 将来のCSV/JSON比較契約

`compare_pile.py` は将来の杭節点変位・断面力・周面摩擦/地盤反力に共通のlong形式を使う。JSONは `{metadata, rows}`、CSVは同じrowsと **明示したJSON metadata sidecar**。各行は厳密に `id,quantity,x_m,y_m,z_m,value`。サンプルは `export_reference.py --verification-only` で生成できる。[軸力CSV](results/axial-fe.csv)、[metadata](results/axial-fe.metadata.json)、[CSV比較結果](results/axial-csv-comparison.json)、[梁JSON比較結果](results/beam-comparison.json)を同梱する。

metadataは `schema,case,units,basis,signs,frame,quantities,expected_ids,provenance` が必須。`case` は `id,model,geometry,materials,loads,boundary,interface,water,reference_state` を必須とし、未知/null/空値は拒否。3Dと1Dは異なるmodelである。`basis` は原点と右手直交3軸を明示。`frame` は工程・荷重率・自重からの増分か総量かを明記。`quantities` は各量の `unit,component,location,definition`、`signs` は各量の正方向、`expected_ids` は量ごとの完全なID集合。`provenance` は元出力のsourceとSHA256。入力ファイル自体のSHAも比較レポートに残す。

両側のmetadataはprovenance以外を完全一致させ、単位系m-kN-kPaと量別単位を確認する。元ファイルの単位・符号・基底を確認してから外部で明示的に変換し、変換根拠/元SHAをprovenanceへ記録する。メタデータを書き換えて違いを隠す運用は許可しない。座標とIDが一致しない場合は自動補間・ID付替え・交集合比較をしない。欠損ID、重複、未知ID/量、NaN/Inf、異なる材料/幾何/荷重/境界/界面/水/基準状態、工程の不一致、符号や基底の不一致を拒否する。異なるメッシュは、別途検証された明示サンプリング/写像ができるまで比較不能。

周面摩擦kPaと地盤反力kN/m、断面力kNと応力kPa、節点変位と面積平均は別量。断面位置・切断面法線・積分/平均方法・共有稜線の扱いもdefinition/location/signsに記載する。3D原出力の実スキーマと断面積分規約は未確認なので、自動抽出器を推測実装しない。比較器は契約に合う将来CSV/JSONを受け取り、量別の絶対最大差/RMS/相対L2を出す。ゼロ参照ノルムは相対差null、絶対差は保持。工学的しきい値は未定義で、差がゼロでも `accepted=null`。正常終了は契約適合と計算完了だけを意味する。

## 試験と再現

新規テストは `checks_pile.py` に置き、既存 `tests/`、pytest収集設定、CI/件数条件は変更しない。標準の `test_*.py` 収集へ追加せず、明示コマンドで実行する。新規44件：閉形式/境界・力とモーメント釣合い・既知指数解・h次数・Hex8線形patch/剛体/エネルギー・元3D体積/荷重・CSV/JSON・不一致拒否。テスト内の丸め検査値や次数検査帯は検算法の実装チェックであり、新しい杭の合否許容値ではない。skipは使わない。

Linux全件件数と実行環境は [検証記録](results/validation.json) とログに記録する。既存Python全件は、現CIと同じ `tests/` と `tools/` の `test_*.py` をunittestで収集する。CI自身は4群に分割しており、今回のLinuxでは2ディレクトリを一括実行する。新規44件と既存全件を混ぜて既存CI件数を増やさない。Windows CIは別の実行結果が確定するまで未確認。タイムアウト/キャンセル/未完了を成功としない。具体的なWindowsコマンド・時間/資源目安・判定は [REPRODUCE.ja.md](REPRODUCE.ja.md) を参照。

Linux新規44件は失敗0・error 0・skip 0。既存PythonはCIと同じ543件を収集し、tests/ 533件（2467.212秒）とtools/ 2件を実行。失敗0・error 2・skip 0、setUpClassエラーにより8件は未実行。既存全件合格とはしない。封印済みNCの平均応力/空間履歴の監査でエラーを検出した。平均応力側は saved mean-stress/complement partition differs from own history、空間側は remeasured JSON differs。両方の詳細を全件ログに保持する。新規杭コードは呼ばれていない。既存コード・テスト・許容値・CI設定を変更せず、この失敗と未実行を保持する。release CLIビルドは成功。前提不足で中断した2回も未完了としてログを保持し、成功に数えない。Windows CIと杭の工学的収束受入は未確認/未達のまま。
