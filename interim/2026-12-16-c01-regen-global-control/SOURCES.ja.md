# 先行事例の確認（実装・再計算に先行）

確認日: 2026-11-14（作業フォルダの日付は依頼された記録識別子）。検索断片と本文確認を区別する。以下は逐語引用ではなく短い要約。文献を本コードの精度証明に使わない。

- Yastrebov, *Numerical Methods in Contact Mechanics* (2013), 著者公開本文 §4.1 の Fig.4.6、§5.2 pp.203–209: 2接線の摩擦円、trial traction、法線圧に依存する半径、trial 方向への radial return を確認。曲線のすべり経路を図示する。今回は円弧・方向転換の比較量（ベクトルすべり、traction、散逸）を選ぶ根拠に使う。§5.1の反復安定性を時間積分精度の証明へ転用しない。https://www.yastrebov.fr/LECTURES/Yastrebov_NMCM_Wiley_ISTE.pdf
- Puso & Laursen, *A mortar segment-to-segment frictional contact method for large deformations*, CMAME 193 (2004) 4891–4913, DOI 10.1016/j.cma.2004.06.001。本文確認は LLNL 著者稿 UCRL-JRNL-200667（2003-10-30）、§1, §6。数値的に客観的なすべりを導入し、開閉・stick/slip の活動制約決定を別の難点として扱う。今回は座標回転整合性と全体開閉対照を別に観測する根拠。mortar の結果を本コードへ直接移植しない。https://digital.library.unt.edu/ark%3A/67531/metadc1410676/m2/1/high_res_d/15013715.pdf
- Abaqus 2026 公式 *FRIC*, DSLIP/SLIP, SED/SFD, DROT: 非回復すべりの増分／累積、弾性エネルギー／増分摩擦散逸、接線基底回転の扱いを区別している。今回は各量の履歴と回転整合性を分離して監査する。DROT の3次元剛体面に関する説明を一般大変形の証明にしない。https://docs.software.vt.edu/abaqusv2026/English/SIMACAESUBRefMap/simasub-c-fric.htm
- Abaqus 2025 公式 *Energy computations in a contact analysis*, energy definitions と Results: 接触ばねの回復可能エネルギー、累積散逸、開閉や法線回転による接触不連続仕事を区別。刻みと力釣合いの精度がエネルギーへ影響する。今回は D と E を別々に評価し、記帳差を全体解の因果説明と同一視しない根拠。https://docs.software.vt.edu/abaqusv2025/English/SIMACAEEXARefMap/simaexa-c-contactenergy.htm
- Abaqus 2026 公式 *Convergence Criteria for Nonlinear Problems*, Severe Discontinuity Iterations: 接触切替に関する判定と平衡残差の両方を扱う。今回は元の残差基準を維持して新規全体平衡を計算し、局所再生で代用しない根拠。文書の刻み制御は今回追加しない。https://docs.software.vt.edu/abaqusv2026/English/SIMACAEANLRefMap/simaanl-c-convergcriteria.htm
- Code_Aster v17 公式 U2.04.04 §1.1: 接触摩擦条件の非微分性と大すべりの幾何学的非線形を区別し、対応付けと拘束問題の解法を分ける。今回は活動状態と局所構成則を分けた観測の補助根拠。Code_Aster の散逸式は今回本文未確認であり引用しない。https://code-aster.org/doc/v17/manuals/man_u/u2/u2.04.04/Introduction.html

確認限界:
- Simo & Laursen (1992), *An augmented Lagrangian treatment of contact problems involving friction*, Computers & Structures 42, 97–116, DOI 10.1016/0045-7949(92)90540-G: 出版社の書誌／要旨は検索で確認、本文取得は403。本文・式番号・精度次数の証拠には使わない。https://www.sciencedirect.com/science/article/pii/004579499290540G
- Wriggers, *Computational Contact Mechanics*, 2nd ed. (2006), DOI 10.1007/978-3-540-32609-0: 出版社の書誌と目次のみ確認。https://link.springer.com/book/10.1007/978-3-540-32609-0
- Laursen, *Computational Contact and Impact Mechanics*, DOI 10.1007/978-3-662-04864-1: 出版社書誌と目次のみ確認（hardcover 2002、同ページ copyright 2003）。https://link.springer.com/book/10.1007/978-3-662-04864-1
- これら2書籍の未確認本文を引用しない。1接線の安定性から2接線の精度・単調性・全体履歴収束を導かない。
