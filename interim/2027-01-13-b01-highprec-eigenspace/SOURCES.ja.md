# 先行事例の調査範囲

## シェルと剛体モード

12-25 の調査者は Lee & Bathe (2004), *Computers & Structures* 82, 945–962 の著者公開本文 §§1–4 を確認し、director、共変ひずみ、tying 補間と物理剛体 6 モードを記録した。Lee, Lee & Bathe (2014), 同誌 138, 12–23 の §§1–2 での等方性・パッチ・零エネルギー検査、Lee ほか (2015), 同誌 153, 148–164 の単一要素・二要素のモード形状の検査も同じ一次資料調査の範囲。今回はその記録を参照した。有限な剛性を持つ低モードについても形状を調べる理由になる。これらの文献の数表や受入値は使用していない。[12-25 の本文確認範囲](../2026-12-25-b01-curved-tri6-cause-diagnosis/primary-source-review.json)

以下は今回本文未取得の書誌。個々の式・数表・試験結果は引用しない。

- R. H. MacNeal & R. L. Harder (1985), “A proposed standard set of problems to test finite element accuracy,” *Finite Elements in Analysis and Design* 1, 3–20. [DOI](https://doi.org/10.1016/0168-874X(85)90003-4)
- K. J. Bathe & E. N. Dvorkin (1986), “A formulation of general shell elements—the use of mixed interpolation of tensorial components,” *International Journal for Numerical Methods in Engineering* 22, 697–722.
- D. J. Allman (1984), “A compatible triangular element including vertex rotations for plane elasticity analysis,” *Computers & Structures* 19. T. J. R. Hughes & F. Brezzi (1989), “On drilling degrees of freedom,” *Computer Methods in Applied Mechanics and Engineering* 72. 本診断では両論文の定式化・係数を採用していない。
- R. D. Cook, D. S. Malkus, M. E. Plesha & R. J. Witt, *Concepts and Applications of Finite Element Analysis*, 4th ed., Wiley. 本文の章・節や固有値試験の数値は未確認。

## 固有空間と摂動

- C. Davis & W. M. Kahan (1970), “The Rotation of Eigenvectors by a Perturbation. III,” *SIAM Journal on Numerical Analysis* 7, 1–46. 本文未取得、書誌のみ。[DOI](https://doi.org/10.1137/0707001)
- Å. Björck & G. H. Golub (1973), “Numerical methods for computing angles between linear subspaces,” *Mathematics of Computation* 27, 579–594. 本文未取得、書誌のみ。[DOI](https://doi.org/10.1090/S0025-5718-1973-0348991-3)
- G. H. Golub & C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press (2013). 本文未取得、書誌のみ。

今回直接読んだ一次ソフトウェア資料は、リポジトリ指定 SciPy 1.18.1 の scipy.linalg.subspace_angles の実装と docstring。直交基底から投影外成分を作り、小さい角をその特異値の arcsin で評価する箇所を確認した。同実装の引用先は Knyazev & Argentati (2002)。本件の角度は独立に保存両基底を直交化し、固定順の投影と高精度 SVD で計算した。mpmath 1.3.0 の eigsy は実対称固有値を昇順で返し、固有ベクトルを列に持つ。実装はインストール物で確認できる。[SciPy の関数](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.subspace_angles.html) ／ [mpmath 固有値実装](https://github.com/mpmath/mpmath/blob/1.3.0/mpmath/matrices/eigen_symmetric.py)

出版社・著者サイトへの今回の取得要求はプロキシで HTTP 403。取得候補 URL は通信記録に限り、未確認 URL と特定の一次本文を対応づけてはいない。[取得記録](retrieval.json)
