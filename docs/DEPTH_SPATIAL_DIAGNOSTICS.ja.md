# 乾燥MC掘削の奥行き差と判定不能の追加診断

PR #13 の最終コミット `13c9c717f58dfa9ad5639487be9d6448e382e6cf` から追加した。main `024835cea9a7c12ad12deda10203f4f45f1fd285` を含む。生産CLI、材料則、数値ゲート、既存の参照値・報告・失敗記録は変更しない。一般3D空間収束、正確なFs、P6は未達である。[元の面内系列](PLASTIC_SPATIAL_CONVERGENCE.ja.md)と[今回の実測記録](../interim/2026-11-07-depth-spatial-diagnostics/README.ja.md)を区別する。

## 奥行きだけを変える別系列

`tools/plastic_spatial_convergence.py --axis depth` は、単位立方体の面内 x/z 分割を各2に固定し、奥行き y を2/4/8分割する。48/96/192 Tet4であり、要素数の立方根からhを作らない。押出しは全域uy=0、有限幅は両端面uy=0。撤去域は x<0.5、z>0.5、有限幅ではさらに y<0.5 に固定する。奥行き細分化で撤去幅が縮む系列にしない。元体積1 m³、撤去後は押出し0.75／有限幅0.875 m³を保持する。

同じ初期化・乾燥非関連MC・施工増分1/16・材料積分幅1e−4を用い、元の力残差1e−8、全体力誤差1e−8、正規化降伏1e−7、工学ひずみ0.05の監査を全段で通過した。変位と応力は従来の交差体積積分で比較する。

|奥行き系列、掘削後|2→4層のRMS差|4→8層のRMS差|条件付き次数|
|---|---:|---:|---:|
|押出し・全変位 mm|0.00797918|0.00448014|0.8327|
|押出し・掘削増分変位 mm|0.00476201|0.00251596|0.9205|
|押出し・全応力 kPa|0.611355|0.346802|0.8179|
|有限幅・全変位 mm|0.01320978|0.00691437|0.9339|
|有限幅・掘削増分変位 mm|0.01855452|0.00737870|1.3303|
|有限幅・全応力 kPa|1.314481|0.953688|判定不能|

有限幅応力の差の比は約0.7255で、従来の0.7条件を満たさない。どの3段階も、固定した面内粗分割に対する奥行き差であり、一般収束・漸近域・連続体解の精度を認定しない。

施工増分1/32・材料積分幅5e−5の同一メッシュ対照を全3段で保存した。最細の全変位RMS差は押出し4.46834e−6／有限幅8.00377e−6 mm、応力差は0.000426011／0.000966976 kPa。空間差より小さいが、全ての積分・非線形誤差が消えた証明ではない。

## 奥行き2層固定の面内系列

押出しの奥行き2層固定で、面内2/4/8の48/192/768 Tet4も新しい出力先へ保存した。全変位RMS差は0.0808553→0.0582971 mm、掘削増分は0.0794616→0.0581268 mm、応力は2.839444→2.004075 kPa。いずれも元の0.7条件に届かず次数は判定不能である。元の1層24/96/384系列の0.081069→0.055366 mm／条件付き次数0.5501や、有限幅系列の既存値を置き換えない。これは別の固定奥行きでの観測であり、異なる系列の差を1本のRichardson列へ混ぜない。

## 押出し応力と掘削モニターの追加診断

`tools/diagnose_plastic_spatial.py` は元の1層面内系列と半減対照を読み、同じ入力SHAによる小規模再測定の全場差を照合する。元の応力差2.787761→1.982204 kPa、比0.7110377722と判定不能はそのままである。応力差のテンソルFrobeniusノルムを成分・静水圧/偏差・塑性履歴の有無で分解した。

|成分|2→4の重み付きRMS kPa|4→8の重み付きRMS kPa|差の比|
|---|---:|---:|---:|
|xx|1.193229|0.914861|0.766710|
|yy|0.756572|0.549328|0.726075|
|zz|2.128856|1.447218|0.679810|
|xy|0.095458|0.057455|0.601891|
|yz|0.268715|0.129081|0.480365|
|xz|1.077993|0.822190|0.762704|

せん断は√2を含む表示で、6成分の二乗和が全応力RMSの二乗になる。静水圧成分と偏差成分の二乗積分の和、塑性履歴分類ごとの二乗積分の和も元の全テンソル量に一致する。分類のplasticは累積相当塑性ひずみ>0という意味で、現在の塑性負荷・接線安定性を意味しない。どれかの成分が減少していても、全応力の元のゲートを置き換えない。

共通位置(0.5,0,0.5)の掘削増分は元の0.140613→0.110561→0.114497 mmを保持する。符号付き差は−0.03005211→+0.00393552 mm。全変位の差−0.04672035→−0.00728537 mmから、初期平衡の差−0.01666824→−0.01122090 mmを引くと、この符号反転が現れる。全変位が単調に近づくことは掘削増分の外挿を保証しない。

半減対照でも0.140610→0.110560→0.114475 mmと振動し、外挿は判定不能である。各段の対照差の絶対値和は、隣接差に対し0.000004888／0.000023059 mmであり、観測した符号反転より小さい。この一回の対照で原因を確定せず、非一意性・局所化・さらに細かい系列での挙動は残件とする。判定不能や不収束を合格、破壊側上端、Fsへ変換しない。

## MCCの既存記録を分離して固定

`tools/summarize_mcc_spatial_records.py` は7つの公開済みJSONの経路とSHA-256を固定し、元の全工程の3種類の変位差を転記する。solverを呼ばず、交差領域積分も再計算しない。欠落・改変はエラーとし、大規模再求解へフォールバックしない。元のstale labelと訂正記録も保持する。

|MCC比較|掘削工程増分の相対L2差|
|---|---:|
|274→945|18.1664%|
|945→1324|12.7877%|
|1324→2262|9.9690%|
|壁固定2262→2865|6.368386%|
|壁固定2865→5144|3.387772%|
|壁固定5144→8061|1.463403911%|
|壁固定8061→14909|1.527772593%|

最後の2対のRMS差も0.0140300→0.0147132 mmと非減少。壁分割の変化がある初期3対と固定壁比較を区別し、今回の乾燥MC面内・奥行き系列へ混ぜない。要素数から次数を推定しない。これらは既存報告の再読込であり、数GBの生結果の新しい独立監査ではない。

既存結果の評価手順は、(1)下の要約コマンドで公開報告のSHAを検査、(2)`from_seed`・`from_construction_origin`・`stage_increment`を全工程で区別、(3)物理初期値/荷重/壁/境界と実行版の元の監査・結合先を確認、(4)撤去後216 m³の共通活動領域と要素別被覆、二乗差・参照積分・絶対RMS・相対L2を一緒に参照、(5)未達を保持、の順に固定する。元報告の結合先は各sourceの`sha256`等に残っている。生結果がローカルに無い場合、この記録以上の再評価はユーザーPCで元の凍結入力・ソルバー・監査一式をそろえて行う。生GP履歴の平均化・移送や新ソルバーでの上書きで代用しない。

## 再現コマンドとユーザーPC向け細分化

固定依存環境でreleaseをビルドし、必ず未使用の出力先を指定する。Windowsの実行ファイルは自動で`.exe`を選択する。

```text
cargo build --locked --release
python tools/plastic_spatial_convergence.py runs/depth-plane-v1 --axis depth
python tools/plastic_spatial_convergence.py runs/depth-finite-v1 --axis depth --finite-width
python tools/plastic_spatial_convergence.py runs/plane-ny2-v1 --depth-divisions 2
python tools/plastic_spatial_convergence.py runs/depth-plane-half-v1 --axis depth --increment 0.03125 --material-substep 0.00005 --control-baseline runs/depth-plane-v1
python tools/plastic_spatial_convergence.py runs/depth-finite-half-v1 --axis depth --finite-width --increment 0.03125 --material-substep 0.00005 --control-baseline runs/depth-finite-v1
python tools/plastic_spatial_convergence.py runs/plane-diagnostic-v1
python tools/diagnose_plastic_spatial.py runs/plane-diagnostic-v1/report.json interim/2026-11-07-plastic-spatial-convergence/plane.json interim/2026-11-07-plastic-spatial-convergence/plane-half.json runs/plane-additional-diagnostic-v1.json
python tools/summarize_mcc_spatial_records.py runs/mcc-record-inventory-v1.json
```

以下の4/8/16と8/16/32は**ユーザーPC専用で、このクラウド環境では実行していない**。奥行き1層押出しの32分割は6144要素、有限幅2層は12288要素。求解と交差領域積分の時間・メモリーに注意し、最後の3段の条件付き次数の持続・応力・モニターを評価する。`--timeout`は求解ごとの壁時計上限であり、交差積分の上限ではない。超過や不収束は未完であり、Fsの上端にしない。

```text
python tools/plastic_spatial_convergence.py runs/pc-plane-4-8-16 --divisions 4 8 16 --timeout 3600
python tools/plastic_spatial_convergence.py runs/pc-plane-8-16-32 --divisions 8 16 32 --timeout 3600
python tools/plastic_spatial_convergence.py runs/pc-finite-4-8-16 --finite-width --divisions 4 8 16 --timeout 3600
python tools/plastic_spatial_convergence.py runs/pc-finite-8-16-32 --finite-width --divisions 8 16 32 --timeout 3600
```

Griffithsの1200/4800/19200 Tet10は既存の入力準備・検査記録を保持する。本タスクでは実求解、661節点の5%比較、正確なFs、実断面SRM、Native UI、構造・水・JOINT、別PC/GPU検証、完成配布は行わない。一般3D空間収束とP6の未達を継続する。
