# 地盤応力の初期化と厚さ方向の積分

## 520cecf から継承した本文確認

Abaqus 2025 の SIGINI はせん断を含む初期応力成分の指定と、静的工程または地盤応力工程で外力との平衡を確認する手順を説明する。Initial Condition Types は土の有効応力の指定と初期平衡、Reactivation は再活性化時の配置基準を扱う。この3件の本文取得・確認は前回に帰属する。今回の比較では本コードの成分順 xx, yy, zz, xy, yz, zx を用いた。[520cecf の一次資料記録](../2027-01-07-a03-explicit-shear-initial-stress/primary-source-review.json)・[SIGINI](https://docs.software.vt.edu/abaqusv2025/English/SIMACAESUBRefMap/simasub-c-sigini.htm)・[Initial Condition Types](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEPRCRefMap/simaprc-c-initialcond-t-SpecifyingTheTypeOfInitialConditionBeingDefined-sma-topic1.htm)・[Reactivation](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEANLRefMap/simaanl-c-elemremovrepl.htm)

## 今回確認した一次資料

- GEOS 開発元の初期化例は、静水圧と重力を与えた準静的な連成平衡を解き、初期化後に変位をゼロに戻す。別の表入力例は、深さに沿って異方的な水平応力を指定した場合、とくに不均質なモデルでは力学的平衡が保証されないと説明する。この2例の応力勾配から本件の明示 τ の参照値は得られない。本件では配置後の保持力と支持拘束を射影した自由残差を別々に測る。[GEOS 重力初期化](https://github.com/GEOS-DEV/GEOS/blob/dee519d6e798ce0c8415e539e9430b625a8aab12/src/coreComponents/physicsSolvers/multiphysics/docs/gravityInducedStressInitialization/Example.rst)・[表入力による初期化](https://github.com/GEOS-DEV/GEOS/blob/dee519d6e798ce0c8415e539e9430b625a8aab12/src/coreComponents/physicsSolvers/multiphysics/docs/userTableStressInitialization/Example.rst)
- OpenSees の LayeredShellFiberSection 実装は各層の厚さから中心位置と重みを作り、層中心の応力に層厚を掛けて断面力を加算する。例えばこの中点則を層中心からの座標 z の二乗に適用すると、厚さ t の区間で中点値による積分は 0、直接積分は t³/12 となる（ここで計算した例）。厚さと評価点の関係を示す断面積分の例で、本件の三次元 Tet10 の誤差評価へ転用しない。[開発元実装](https://github.com/OpenSees/OpenSees/blob/5b9b70511e04e8f1b65710edb98fdffd8ef3add0/SRC/material/section/LayeredShellFiberSection.cpp)
- SciPy の roots_legendre の開発元本文は、n 点 Gauss–Legendre 則の多項式厳密次数 2n−1 を標準区間・重み 1 について記している。今回測る GL5−GL7 は厚さを変えた有限要素の体積再積分・再求解の差である。そこから連続体解の誤差や一般の薄層の空間収束次数は決定できない。[SciPy v1.15.2](https://github.com/scipy/scipy/blob/v1.15.2/scipy/special/_orthogonal.py#L2319)

Abaqus/Standard “Geostatic Stress State” 2024・2025 年版と 2025 年版の shell section behavior は HTTPS CONNECT 403 で本文未取得、書誌のみ。今回本文を確認した地盤例には明示 τxz と同じ初期値での平衡解はない。外部参照値は導入しない。[Geostatic Stress State](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEANLRefMap/simaanl-c-geostatstress.htm)・[Shell Section Behavior](https://docs.software.vt.edu/abaqusv2025/English/SIMACAEELMRefMap/simaelm-c-shellsectionbehavior.htm)
