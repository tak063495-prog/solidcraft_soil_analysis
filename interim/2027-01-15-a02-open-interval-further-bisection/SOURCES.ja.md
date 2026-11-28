# 事象位置決め・増分制御の一次資料

2026-11-20 JST に下記の公開本文・ソースを参照した。2027-01-15 は依頼されたフォルダ名の日付。外部の数値を検証値に採用していない。

- SUNDIALS CVODE の “Rootfinding” は利用者が指定した関数 g(t,y(t)) の符号変化を検出し、区間内で修正割線法を使って根の位置を求める。符号を変えない偶数重根を見逃し得る。本診断が保持するのは、入力を細分し全体平衡を解き直した各受理終点とその直前の終点である。この観測には同一の連続軌道上で g を評価する操作がない。[開発元マニュアル本文](https://raw.githubusercontent.com/LLNL/sundials/main/doc/cvode/guide/source/Mathematics.rst)
- deal.II の step-42 の接触例では、法線の乗数と変位・初期 gap の条件から接触自由度を active set に加え、次の Newton 更新に与える拘束を組む。GetFEM の接触摩擦の解説は法線の片側条件と、接線力の摩擦円盤への射影を記述する。これらが説明するのは現在の変位・力による状態判定であり、本例で受理点の間を通る時刻を与える式は引用していない。[deal.II 実装](https://raw.githubusercontent.com/dealii/dealii/master/examples/step-42/step-42.cc)・[GetFEM マニュアル本文](https://raw.githubusercontent.com/getfem-doc/getfem/master/doc/sphinx/source/userdoc/model_contact_friction.rst)
- ANSYS PyMAPDL の AUTOTS の Notes は、自動増分 ON で刻み予測と二分を併用すると説明する。NSUBST の初期サブステップ数は接触要素の有無にも依存する。v242 の Contact Technology Guide は直接取得が HTTP CONNECT 403 で書誌のみ。別途 v252 の公式本文 §3.9.21.1–2 を確認した。KEYOPT(7)=1 では反復中の過大貫入、適応 small-sliding での過大すべり距離、大きな接触状態変化が二分の条件。次 substep の増分予測は KEYOPT(7)=2,3 で別に記述されている。[Contact Technology Guide v252](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_ctec/Hlp_ctec_realkey.html) [ANSYS コマンド本文](https://raw.githubusercontent.com/ansys/pymapdl/main/src/ansys/mapdl/core/_commands/solution/load_step_options.py)
- Abaqus/Standard の FRIC の PNEWDT と DSLIP の本文確認は 01-05 の一次資料取得記録を継承する。局所から再試行を要求する操作と、収束増分で不可逆すべりを累積する操作を区別する。今回、接触診断・収束制御の公式ページの直接取得は HTTP CONNECT 403。この 2024 版の2ページは書誌のみ。別途 2026 版 “Common Difficulties” の貫入距離による増分制御の本文を確認した。適用対象の接触では未収束反復中の貫入距離が限界を超すと増分を破棄し、より小さい増分で再試行する。有限すべりの面対面接触など、同じ規定が適用されない場合も明記されている。[Common Difficulties](https://docs.software.vt.edu/abaqusv2026/English/SIMACAEITNRefMap/simaitn-c-contacttrouble.htm) [FRIC](https://docs.software.vt.edu/abaqusv2026/English/SIMACAESUBRefMap/simasub-c-fric.htm)・[01-05 の取得記録](../2027-01-05-a02-n8-contact-cap-open-interval/primary-source-review.json)・[Contact Diagnostics](https://docs.software.vt.edu/abaqusv2024/English/SIMACAEITNRefMap/simaitn-c-contactdiagnostics.htm)
- Hairer, Nørsett & Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed., Springer. 出版社サイトへの直接取得は失敗。その後、原著の第 II.6 節 pp. 195–196 の本文を公開複本で確認した。Algorithm 6.4 は事象関数の符号変化を探し、dense output を代入した方程式を二分法や Newton 法で解く。1刻み中に2度符号が変わる場合は端点だけで検出できない。“Discontinuous Equations” の式 (6.23) は、切替関数の正負で異なる右辺を与える。ここでの入力格子二分は各分割で接触履歴自体を更新するので、この根探索の収束をそのまま適用できない。[目次と頁の照合](https://toc.library.ethz.ch/objects/pdf/e01_978-3-642-05163-0_01.pdf)・[著者の書誌](https://www.unige.ch/~hairer/books.html)・[DOI](https://doi.org/10.1007/978-3-540-78862-1)

今回取得した本文の SHA-256（取得した順に依存しない）。

- CVODE Mathematics.rst　`2cbd13b088bfb8627477934f9c8c175e6f0731b68cbe2e7f0e9bb845c93f1e1c`
- deal.II step-42.cc　`c8e1bdbc7161592ff7b50943792e73c6aa04cc1d2bb237649ac433223701012f`
- GetFEM model_contact_friction.rst　`2bf1a317a96f679d5cef7ee888df013f1ba040979955fa07bf6223dde29b22c9`
- ANSYS load_step_options.py　`cedb3a0e881d4caa6532cb72f5c24e09cb4c8429d2555003a2bc962753f68453`
