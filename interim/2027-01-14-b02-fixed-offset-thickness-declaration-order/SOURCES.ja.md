# 板厚・腕長と演算順

OpenSees の RigidBeam は保持節点と拘束節点の座標差から小回転に対する並進拘束行列を作る。ElasticMembranePlateSection は厚さ h を引数として受け、せん断剛性に h、曲げ剛性に h³ を使う。両ソースで幾何学的な腕と断面厚を別々に確認した。今回も腕長 e を m で固定し、板厚変更に伴う剛性・自重の変化を含めて観測する。[RigidBeam.cpp](https://github.com/OpenSees/OpenSees/blob/5b9b70511e04e8f1b65710edb98fdffd8ef3add0/SRC/domain/constraints/RigidBeam.cpp#L124-L155) ／ [ElasticMembranePlateSection.cpp](https://github.com/OpenSees/OpenSees/blob/5b9b70511e04e8f1b65710edb98fdffd8ef3add0/SRC/material/section/ElasticMembranePlateSection.cpp#L195-L228)

Abaqus の参照面オフセットは板厚の符号付き比で指定する。同じ物理的な腕長を保つ場合の比は e/t に応じて変わる。Abaqus と本リポジトリ間の法線・符号規約の一致は仮定しない。これは 12-24 に公式本文を確認した記録の引用。今回の公式サイトへの接続は 403 だった。[12-24 の一次資料台帳](../2026-12-24-b02-offset-composite-cause-diagnosis/primary-source-review.json) ／ [Abaqus, Defining the Initial Geometry of Conventional Shell Elements](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEELMRefMap/simaelm-c-shellgeometry.htm)

Ansys の PyDyna が公開する SECTION_SHELL の NLOC 説明では、節点面から参照面までのシェル法線方向距離を −0.5 × NLOC × 平均板厚とする。したがって NLOC 一定で厚さを変える操作では距離も変わる。確認した本文は Ansys のカード API 説明であり、LS-DYNA マニュアル本体は未取得。符号や距離を本モデルに転用していない。[Ansys PyDyna, SectionShell.nloc](https://github.com/ansys/pydyna/blob/553c61c429c61fc759492012a5de1b364666dbcb/src/ansys/dyna/core/keywords/keyword_classes/auto/section/section_shell.py#L344-L350)

Open MPI の MPI_Reduce マニュアルは、評価順の変更で浮動小数点加算の結果が変わり得ると説明している。今回調べるのは宣言配列の順序。節点 ID、要素内接続、物理的な自由度と成分の対応を固定し、追加の自由度再番号付け試験は行わない。各実装の内部加算・連立求解・停止判定への寄与を、この比較だけから分配することはできない。[Open MPI, MPI_Reduce](https://github.com/open-mpi/ompi/blob/363df5354f9481830a76cba24000e22a1f8562d5/docs/man-openmpi/man3/MPI_Reduce.3.rst#L109-L117)

本文未取得の書誌・項目名は Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM (2002)、Nastran のシェル ZOFFS、LS-DYNA Keyword User’s Manual の *SECTION_SHELL（NLOC）。これらの式・符号・適用範囲や文献の参照数値は使用していない。[Higham 書誌](https://doi.org/10.1137/1.9780898718027)
