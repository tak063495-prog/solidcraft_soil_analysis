# T9 出典と確認範囲

2026-11-22 に照合した公開資料。番号は [仕様調査](README.ja.md) の各判断に対応する。製品資料はその手続きの意味、ベンチマークは比較設計の根拠として使った。

|番号|書誌・確認箇所|確認した内容と適用範囲|
|---|---|---|
|1|Seequent / Bentley Systems, *PLAXIS 3D 2025.1 Reference Manual*, §§7.6.1.2, 7.9.1.4。[公式 PDF](https://files.seequent.com/PLAXIS/Manuals/PLAXIS_3D/English/PLAXIS_3D_2_Reference%20Manual.pdf)|段階施工の目標、工程荷重の内挿、未完の三つの理由。2D 版でも崩壊表示に用いる連続除荷と現剛性指標 CSP を確認した。[2D Reference Manual 2025.1](https://files.seequent.com/PLAXIS/Manuals/PLAXIS_2D/English/PLAXIS_2D_2_Reference%20Manual.pdf)|
|2|Dassault Systèmes, *Abaqus 2026 Analysis Guide*, “Element and Contact Pair Removal and Reactivation” および “Unstable Collapse and Postbuckling Analysis”。[要素撤去](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEANLRefMap/simaanl-c-elemremovrepl.htm)・[Riks 法](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEANLRefMap/simaanl-c-postbuckling.htm)|境界節点力の漸減、static Riks ステップの撤去制限、比例荷重の経路追跡。|
|3|code_aster, *DocAster*, U2.04.06 §4（掘削モデル・節点反力）、R5.03.80 §2（荷重制御）。[掘削](https://codeaster.gitlab.io/doc/docaster/manuals/man_u/u2/u2.04.06/Avant_de_s_attaquer_au_ficher_de_commande_Code_Aster_.html)・[追跡法](https://codeaster.gitlab.io/doc/docaster/manuals/man_r/r5/r5.03.80/Principe_des_m_thodes_de_pilotage_du_chargement.html)|トンネル境界の反力と解放、変位条件に従った荷重係数の決定。|
|4|DIANA FEA, *DIANA 10.10 User's Manuals*, “Phased Analysis”。[マニュアル](https://manuals.dianafea.com/d1010/en/1353862-1354346-phased-analysis.html)|要素・支持条件の変更、初期状態の継承、工程ごとの増分荷重。|
|5|Itasca, *FLAC3D 7.0 Documentation*, “zone relax command”。[コマンド仕様](https://docs.itascacg.com/flac3d700/flac3d/zone/doc/manual/zone_manual/zone_commands/cmd_zone.relax.html?highlight=zone+relax)|掘削区画の漸減、力の比による既定のサーボ制御、急撤去に伴う疑似慣性の抑制。|
|6|C. McGann & P. Arduino, University of Washington, “Excavation Supported by Cantilevered Sheet Pile Wall”, *OpenSees Examples*, Model Description / Boundary Conditions。[解説・入力ファイル](https://opensees.berkeley.edu/wiki/index.php/Excavation_Supported_by_Cantilevered_Sheet_Pile_Wall)|側方拘束を持つ二次元例。層を外し、平衡を回復してから次層へ進む。|
|7|M. A. Crisfield (1981), “A fast incremental/iterative solution procedure that handles ‘snap-through’”, *Computers & Structures* 13(1–3), 55–62. DOI 10.1016/0045-7949(81)90108-5。[出版社抄録](https://www.sciencedirect.com/science/article/abs/pii/0045794981901085)|荷重を未知数とし、増分の長さの制約を加えて極限点を通過する原理。|
|8|R. de Borst, M. A. Crisfield, J. J. C. Remmers & C. V. Verhoosel (2012), *Non-Linear Finite Element Analysis of Solids and Structures*, Wiley, ch. 4, “Solution Techniques in Quasi-Static Analysis”。[出版社の章概要](https://onlinelibrary.wiley.com/doi/10.1002/9781118375938.ch4)|準静的解析での経路追跡・弧長法の位置づけ。|
|9|D. M. Potts & L. Zdravković, *Finite Element Analysis in Geotechnical Engineering*, Theory and Application, Thomas Telford（2巻セット、2001）。[出版社紹介](https://bookstore.emerald.com/media/flyers/9780727728128.pdf)。D. M. Potts & A. B. Fourie (1984), “The behaviour of a propped retaining wall: results of a numerical experiment”, *Géotechnique* 34(3), 383–404。[論文抄録](https://doi.org/10.1680/geot.1984.34.3.383)|書籍はモデルの仮定と実適用の位置づけを確認。初期応力と施工方法が壁応答に与える影響の具体例は後者の論文による。|
|10|Dassault Systèmes, *Abaqus 2025 Benchmarks Guide*, “3DNLG-7: Elastic large deflection response of a hinged spherical shell under pressure loading”, Reference solution / Results and discussion. 出典は NAFEMS R0024, test 3DNLG-7。[比較表と参照解の由来](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-3dnlg7.htm)|圧力と中央変位を二つの極限点で比較。掲載参照値の生成元に関する独立性の注意も明記されている。|
|11|PETSc Development Team, *PETSc 3.24.6 Users Manual*, “KSP: Linear System Solvers”, Convergence Tests。[版固定のマニュアル](https://web.cels.anl.gov/projects/petsc/vault/petsc-3.24.6/docs/manual/ksp.html)|相対・絶対許容値、右辺の大きさ、前処理の有無とノルムの選択。|
|12|R. Barrett et al. (1994), *Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods*, SIAM, §4.2 “Stopping Criteria”。[Netlib 配布版](https://www.netlib.org/templates/templates.pdf)|残差、後退誤差、逆行列ノルムによる前進誤差の評価。|
|13|LAPACK, *DGESVX*, double-precision general linear system expert driver, Purpose / Algorithm。[Netlib 公式ソースと仕様](https://www.netlib.org/lapack/explore-html/da/d21/dgesvx_8f_source.html)|平衡化、ピボット付き LU、逆条件数推定、反復改良、前進・後退誤差の出力。|
|14|D. V. Griffiths & N. Koutsabeloulis (1985), “Finite element analysis of vertical excavations”, *Computers and Geotechnics* 1(3), 221–235。[出版社抄録](https://www.sciencedirect.com/science/article/pii/0266352X85900254)|平面ひずみ・軸対称掘削の限界高さを遠心模型試験・極限解析と照合し、高さ／幅比に伴う機構の違いを評価。|

## プロジェクト固有の一次記録

- 案1の共通 λ による局所経路と左右弱モード。[b069ec8 の診断](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/b069ec84da432f677c33fab5c984b6e680835bc4/interim/2026-11-17-nc-fold-stability)
- 荷重源、全活動積分点、非一様場の監査。[11月13日・物理限界の検討](../2026-11-13-nc-lambda-physical-limit/README.ja.md)
- 固定強度の離散上下側値、条件掃引、折返し対分岐の分類 D。[11月13日・崩壊分類](../2026-11-13-nc-lambda-collapse-classification/README.ja.md)
- PC 実行結果（ユーザーPC、未コミット）。2026-11-22、実行時先端 e6181bae、提供元の保存場所は `/workspace/pc_runs/fine_20261122/nc96/`。本調査の PC 数値はユーザーが提示した記録の転記であり、原ログの再集計は行っていない。
