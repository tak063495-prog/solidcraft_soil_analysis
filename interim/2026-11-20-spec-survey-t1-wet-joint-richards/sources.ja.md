# 湿潤 JOINT と Richards 連成の出典

確認日 2026-11-20。README の角括弧内の番号は下表に対応する。「本文」は該当する節・式を閲覧、「要旨」は出版社の要旨までを確認したもの。製品の年が判明しない場合は版を識別子とした。シェルから出版社等への直接取得は HTTP 403。Web 検索・閲覧で以下を照合した。製品名と本コードのモード名の対応は本調査の提案である。

|番号|著者・年、題名、誌名またはマニュアル・該当箇所|確認状況と引用範囲|
|---|---|---|
|1|J. M. Segura, I. Carol（2008年オンライン刊行）, “Coupled HM analysis using zero-thickness interface elements with double nodes. Part I: Theoretical model”, *International Journal for Numerical and Analytical Methods in Geomechanics*, Abstract。[DOI](https://doi.org/10.1002/nag.735)|出版社要旨。変位と圧力の節点未知量、開口依存の透水性、有効応力と相対変位、逐次と一体 Newton、非対称 Jacobian を確認。RT0 面写像の先例には使わない。|
|2|J. M. Segura, I. Carol（2008年オンライン刊行）, 同 Part II: “Verification and application”, 同誌, Abstract。[DOI](https://doi.org/10.1002/nag.730)|出版社要旨。鉛直不連続面を加えた圧密で staggered と fully coupled を比較。貯留層内の水平接合の例もある。検証の形状・目的だけを参照する。|
|3|D. Garolera, J. M. Segura, I. Carol, M. R. Lakshmikantha, J. Alvarellos（2017）, “Coupled H-M fracture interaction using FEM with zero-thickness interface elements”, *COMPLAS XIV*, §2、式(7)〜(12)。[著者所属大学の会議論文](https://upcommons.upc.edu/bitstreams/6b5922db-6255-469f-8930-f6928a9c1c5d/download)|本文。両面節点圧を中面に補間する平均、上下の圧力差、力学弱形式を確認。平均と差に別の役割がある。|
|4|Dassault Systèmes SIMULIA（2025版）, *Abaqus Elements Guide*, “Defining the Constitutive Response of Fluid within the Cohesive Element Gap”, Tangential Flow within the Gap / Normal Flow across Gap Surfaces。[マニュアル](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEELMRefMap/simaelm-c-cohesivefluidbehavior.htm)|本文。Newton 流体の開口三乗係数、中面圧と各面圧の差による両側への漏出、開口変化を含む連続式。充水した cohesive 要素を対象とする。|
|5|Bentley Systems, *PLAXIS 3D Reference Manual*, CONNECT Edition V22.02, “Modelling loads and structures”, pp.127–128 表3、および “Groundwater flow through interfaces” p.259参照。[公式版](https://bentleysystems.service-now.com/sys_attachment.do?sys_id=fadeca1f9751661081d373b0f053af13)|表3本文。相を通じた Active in flow の設定、横断抵抗と面内 drainage conductivity の区別。開閉から透水性が自動的に決まるとの記載を根拠にしていない。|
|6|DIANA FEA BV, *DIANA Theory Manual* 10.4, “Detailed Groundwater Flow”, §59.2.6 Resistance Layers。[理論書](https://manuals.dianafea.com/d104/Theory/Theoryse380.html)|本文。横断する抵抗層は水理伝導度を物理厚さで割った係数、単位の例は s⁻¹。圧力水頭への依存を許す。|
|7|DIANA FEA BV（2022更新、10.8版）, *DIANA Theory Manual*, “Pore pressure and Hydraulic Heads”, Direct Pore Pressure / Fluid Pore Pressure / Hydraulic Head, 式(1)〜(8)。[理論書](https://manuals.dianafea.com/d108/en/1219784-1219845-pore-pressure-and-hydraulic-heads.html)|本文。先行流れ解析からの分布荷重、界面法線荷重と飽和度の扱いを確認。式(8)の正圧への制限を本コードの土の負圧に適用しない。|
|8|Itasca, *FLAC3D 7.0*, Interfaces, “Formulation” および “Modeling Guidelines”。[構成則](https://docs.itascacg.com/flac3d700/flac3d/docproject/source/theory/interfaces/formulation.html)・[適用条件](https://docs.itascacg.com/flac3d700/flac3d/docproject/source/theory/interfaces/modelingguidelines.html)|本文。target 面からの圧力補間と有効法線力に基づくすべり。ここでの水圧は面を押し開く力を与えない。|
|9|Itasca, *3DEC 7.0*, “Joint Fluid Flow”, Hydraulic Aperture / Running the Coupled Solution。[理論書](https://docs.itascacg.com/3dec700/3dec/docproject/source/theory/fluid/fluid_joint.html)|本文。接触変位からの水理開口、圧力と力学の更新、機械・流体ステップの組み合わせ。FLAC3D の interface とは別の実装。|
|10|EDF, *Code_Aster*, R7.01.25, “Lois de comportement des joints des barrages”, §5, §6。[規定圧](https://codeaster.gitlab.io/doc/docaster/manuals/man_r/r7/r7.01.25/Prise_en_compte_de_la_pression_hydrostatique_sans_couplage.html)・[HM 連成](https://codeaster.gitlab.io/doc/docaster/manuals/man_r/r7/r7.01.25/Formulation_th_orique_du_couplage_hydrom_canique.html)|本文。規定圧の法線応答と XXX_JOINT_HYME の力学・流れ、開口依存の流量、最小開口を区別する。不飽和過渡や世代管理への適用は未確認。|
|11|P. A. Witherspoon, J. S. Y. Wang, K. Iwai, J. E. Gale（1980）, “Validity of Cubic Law for fluid flow in a deformable rock fracture”, *Water Resources Research* 16(6), 1016–1024, Abstract。[DOI](https://doi.org/10.1029/WR016i006p01016)|出版社要旨。岩石割れ目の実験と平行平板の三乗則。粗さや不飽和の土・構造接触への一般適用、残留開口値は未確認。|
|12|M. A. Celia, E. T. Bouloutas, R. L. Zarba（1990）, “A general mass-conservative numerical solution for the unsaturated flow equation”, *Water Resources Research* 26(7), 1483–1496, Abstract。[DOI](https://doi.org/10.1029/WR026i007p01483)|出版社要旨。Richards 式の保存形と水収支、解の精度を区別する根拠。界面水量の許容値は与えない。|
|13|J. Kim, H. A. Tchelepi, R. Juanes（2011）, “Stability and convergence of sequential methods for coupled flow and geomechanics: Fixed-stress and fixed-strain splits”, *Computer Methods in Applied Mechanics and Engineering* 200, 1591–1606, Abstract・安定性解析。[DOI](https://doi.org/10.1016/j.cma.2010.12.022)・[著者配布版](https://juanesgroup.mit.edu/publications?action=AttachFile&do=get&target=kim_tchelepi_juanes_cmame_2011-fixedstress.pdf)|著者版本文・要旨。逐次法では分割が安定性を左右する。非線形摩擦開閉を含む本系の収束保証は未確認。|
|14|Dassault Systèmes SIMULIA（2025版）, *Abaqus User Subroutines Guide*, “FRIC”, DSLIP / SLIP / SFD。[マニュアル](https://docs.software.vt.edu/abaqusv2025/English/SIMACAESUBRefMap/simasub-c-fric.htm)|本文。不可逆すべり増分が出力判定に用いられ、収束時に累積される。物理的に何も変えない記録フレームで旧ラベルを保持する規約は未確認。|
|15|Dassault Systèmes SIMULIA（2025版）, *Abaqus Analysis Guide*, “Element and Contact Pair Removal and Reactivation”, strain-free reactivation / material state。[マニュアル](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEANLRefMap/simaanl-c-elemremovrepl.htm)|本文。無ひずみ再活性化、多孔質材料の間隙率は初期値、飽和度は除去時の値を保持する例。これは新設湿潤 JOINT の水源・質量移送の規約ではない。|
|16|R. W. Lewis, B. A. Schrefler（1998）, *The Finite Element Method in the Static and Dynamic Deformation and Consolidation of Porous Media*, 2nd ed., Wiley, ISBN 9780471928096。[出版社](https://www.wiley.com/en-gb/the-finite-element-method-in-the-static-and-dynamic-deformation-and-consolidation-of-porous-media-2nd-edition-p-9780471928096)|書誌・出版社による内容紹介を確認。部分飽和多孔質体の連成の背景。JOINT に関する個別節・式の本文は未確認。|

## リポジトリの一次記録

以下は起点 `82409c56f80022d3d161f9d0ef32afa19bfa3cfa` の記録。ディレクトリの日付は保存ラベルである。

|番号|根拠と確認箇所|
|---|---|
|R1|[issue triage の todo](../2026-11-16-issue-triage/todo.ja.md) の A-01/A-04/B-02/C-04 と追補。[A-04 検証契約](../2026-11-29-richards-wet-joint-a04/SPEC.ja.md)・[A-04 受理範囲・不整合の拒否漏れ](../2026-11-29-richards-wet-joint-a04/README.ja.md)。|
|R2|[A-04 S1〜S12 の継続症状と未達理由](../2026-12-23-a04-richards-wet-joint-cause-diagnosis/README.ja.md)、[S1〜S12 追補](../2026-12-31-a04-richards-wet-joint-cloud-followup/README.ja.md)、[第4反復の乾燥対照と仕様判断待ち](../2027-01-17-a04-iteration4-tangent-rhs-correction/README.ja.md)。|
|R3|[C-04 規定水位と単独 PDE](../2026-11-22-wet-joint-transient-water/README.ja.md)、[C-04 16/64 と厳密カテゴリ](../2027-01-03-c04-transient-water-cloud-followup/README.ja.md)。|
|R4|[A-01 組合せ、世代と再開](../2026-11-23-joint-combo-extension/README.ja.md)、[A-01 残件](../2027-01-04-a01-high-precision-nullspace/README.ja.md)、[B-02 受理表](../2026-12-01-tri6-offset-composite-b02/README.ja.md)、[B-02 の作用面・仕事の判断](../2027-01-14-b02-fixed-offset-thickness-declaration-order/README.ja.md)。|
|R5|[P6-01 台帳への追加手順](../2026-12-10-ci-p601-inventory-fix/README.ja.md) の「今後、新規 interim フォルダを追加するとき」。|

OpenSees の乾燥接触要素は公開されているが、ここで問う Richards 面写像・増分ゼロの保持・湿潤 BIRTH の一式について一次資料による採用例は未確認。ZSoil、MIDAS GTS、GeoStudio 固有の同等契約も未確認。これらを判定の根拠に含めていない。
