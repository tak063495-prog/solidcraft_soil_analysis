# 乾燥MC掘削の面内h細分とGriffiths斜面の準備

2026-11-07。一般3D空間収束とP6は未完了のままである。数値基準、許容値、参照値、過去の失敗記録は変更していない。今回の範囲は生産CLIを使う小規模な乾燥MCの掘削と、ユーザーPC向け斜面入力の準備であり、GUI/GPUは対象外である。

## 既存記録から分けること

[空間収束ハーネス](SPATIAL_CONVERGENCE.ja.md)の平滑弾性と一様応力Fs=4/3は、非一様な塑性機構の収束を証明しない。[交差領域の掘削比較](MCC_MESH_OVERLAP.ja.md)では、平均底面隆起が近づいても全場の差が減らない例がある。274→945→1324→2262要素の掘削増分L2差18.1664→12.7877→9.9690％、固定壁の6.368386％、5,144→8,061→14,909要素の非減少は未達のまま残す。これらは今回の比2の面内系列へ混ぜず、要素数の立方根だけからhや次数を作らない。

[平面ひずみ比較](PLANE_STRAIN_2D.ja.md)に保存されたGriffiths斜面の15.582／15.886／15.593／15.672％は、奥行き、施工刻み、材料積分幅を変更した記録であり、面内h系列ではない。元の661節点変位データは本リポジトリに無く、5％比較は今回も完了にしない。

## 小規模ケースと測定

`tools/plastic_spatial_convergence.py` は同じ1 m立方体をTet4で解く。E=10000 kPa、ν=0.25、γ=20 kN/m³、c=2 kPa、φ=20°、ψ=0°、非関連流れ、Rankine面なし。種応力・塑性履歴は0で、自重による初期釣合いから掘削する。側面x変位と底面z変位を固定する。

- 押出し：全節点のy変位を固定し、x<0.5 mかつz>0.5 mを撤去する。奥行きは1層で固定。活動体積0.75 m³。
- 有限幅：y=0/1 m端面のy変位を固定し、内部のy変位は自由。上記領域のうちy<0.5 mだけを撤去する。奥行きは2層で固定。活動体積0.875 m³。この系列は奥行き収束の証明ではない。

x/z分割2/4/8で、要素数は押出し24/96/384、有限幅48/192/768。荷重履歴、材料、境界の物理的位置、撤去領域は同じである。施工増分1/16、最小増分1/4096、材料積分幅1e−4、力残差1e−8、線形許容1e−10を固定する。既存斜面プロトコルの力残差1e−8・全体力1e−8・正規化降伏1e−7・工学ひずみ0.05のゲートを満たした完成Runだけを比較する。小規模ケースの各求解は30秒以内で終了した。

比較は既存の `tet10_overlap.Overlay` の交差多面体と要素別被覆検査を使う。Tet4の一次変位を中点で厳密に二次表現へ載せるのは積分のためだけであり、求解メッシュや結果は変えない。全変位と初期釣合いからの施工増分を分け、同じ活動領域の体積RMSと相対L2を保存する。Tet4の応力は要素ごとに一定なので、交差体積で全6成分の差を厳密に積分できる。応力テンソルのFrobeniusノルムではせん断3成分の自乗に2を掛ける。履歴の移送・応力の平滑化は行わない。

次数は `log2(||q_h−q_h/2|| / ||q_h/2−q_h/4||)`。スカラーは最後の差の符号が同じで、差が0.7倍未満へ減り、次数が既存の[0.5,8]に入る場合だけ参照なしRichardsonを出す。変化が値の1e−8以下の場合は力の求解許容と分離できる次数を付けない。定数に近い平均σzzから次数を作らない。場の差のノルムから出す最細誤差 `Dfine/(2^p−1)` は同じ主要誤差場を仮定する条件付き推定であり、厳密な誤差上界ではない。3段階だけでは漸近域の持続を証明できない。

## 結果

掘削後の体積積分による差。左右は2→4分割と4→8分割である。

|系列|量|粗中差|中細差|条件付き次数|最細誤差推定|
|---|---|---:|---:|---:|---:|
|押出し|全変位RMS mm|0.081069|0.055366|0.5501|0.119264 mm|
|押出し|掘削増分RMS mm|0.079655|0.055181|0.5296|0.124415 mm|
|押出し|全応力RMS kPa|2.787761|1.982204|判定不能|なし|
|有限幅|全変位RMS mm|0.036016|0.016846|1.0963|0.014802 mm|
|有限幅|掘削増分RMS mm|0.029614|0.016135|0.8761|0.019313 mm|
|有限幅|全応力RMS kPa|2.507640|1.672348|0.5845|3.348226 kPa|

押出し応力の差の比は約0.711で、0.7条件を満たさない。減少していることだけを記録し、外挿しない。有限幅の応力推定誤差も大きく、一般収束や実用精度の受入へ置き換えない。

共通の物質位置(0.5,0,0.5) mでのz変位と、活動領域の平均σxxをスカラーとして別に記録した。

|系列・量|2分割|4分割|8分割|次数|連続側の条件付き推定|最細の離散成分|
|---|---:|---:|---:|---:|---:|---:|
|押出し・全uz mm|−0.447557|−0.494277|−0.501563|2.6810|−0.502909 mm|0.001346 mm|
|有限幅・全uz mm|−0.470989|−0.494400|−0.505162|1.1212|−0.514319 mm|0.009157 mm|
|押出し・平均σxx kPa|2.622973|3.000113|3.240187|0.6516|3.660685 kPa|0.420498 kPa|
|有限幅・平均σxx kPa|3.038021|3.156980|3.220623|0.9024|3.293845 kPa|0.073222 kPa|

押出しの同じ位置の掘削増分は0.140613→0.110561→0.114497 mmと振動するため、スカラーの外挿は判定不能である。局所の全変位の高い次数を、掘削増分の場の精度へ読み替えない。掘削後の活動塑性体積は押出し0.250000→0.281250→0.328125 m³、有限幅0.166667→0.218750→0.257813 m³で、弾性だけのパッチ試験ではない。

施工増分1/32・材料積分幅5e−5へ同時に半減した全3段階の対照も保存した。最細の全変位RMS差は押出し1.83796e−5 mm、有限幅3.03084e−6 mm、全応力RMS差はそれぞれ0.00136215／0.000433763 kPaで、上の空間差より小さい。空間差とこの制御変更の差を混ぜず、同じメッシュ上で比較する。全数値は `increment-controls.json` に保存する。この一回の対照だけで全ての積分・非線形誤差が消えたとは扱わない。

非一様な掘削のFsは、独立に釣合いを排除した上端を今回用意していないので探索していない。不収束を破壊側上端にしない。一様応力Fs=4/3の既存結果は変更していない。

## 再実行とユーザーPCでの細分化

固定依存環境とrelease CLIを用いる。出力先は新しいディレクトリを指定する。Windowsでは `.exe` を自動選択する。CIは通常のPython unittestから小規模ケース、最細段の半減対照、既知場と拒否例を実行する。

```text
cargo build --locked --release
python tools/plastic_spatial_convergence.py runs/plastic-plane-ci
python tools/plastic_spatial_convergence.py runs/plastic-finite-ci --finite-width
python tools/plastic_spatial_convergence.py runs/plastic-finite-half --finite-width --increment 0.03125 --material-substep 0.00005 --control-baseline runs/plastic-finite-ci
```

以下はこの環境では実行していない。ユーザーPCで4/8/16の後、8/16/32へ進め、最後の3段で次数が維持されるか確かめる。32分割の要素数は押出し6,144、有限幅12,288。求解と交差領域積分の時間・記憶量は増える。`--timeout` は各求解の壁時計上限で、超過は未完として記録し、その系列からFsや収束合格を作らない。交差領域積分の時間上限ではない。

```text
python tools/plastic_spatial_convergence.py runs/plastic-plane-user --divisions 4 8 16 --timeout 3600
python tools/plastic_spatial_convergence.py runs/plastic-finite-user --finite-width --divisions 4 8 16 --timeout 3600
python tools/plastic_spatial_convergence.py runs/plastic-finite-finer --finite-width --divisions 8 16 32 --timeout 3600
python tools/plastic_spatial_convergence.py runs/plastic-finite-user-half --finite-width --divisions 4 8 16 --increment 0.03125 --material-substep 0.00005 --timeout 3600 --control-baseline runs/plastic-finite-user
```

`report.json` に各段の入力・結果SHA、ソルバーSHA、RMS・相対L2、塑性体積、力・降伏・ひずみ、次数・外挿・拒否理由を保存する。途中の `result.json` は生の履歴を残す。既存出力の上書きはしない。

## Griffiths斜面の面内系列

`tests/data/griffiths_corner_geometry.json` は固定版 `440d246d9e32f5892c79539d4a7b797ad29442ee` の `docs/inconsistent-evidence24/fixtures.json`（SHA-256 `cb706a29d979c0552126baf3107cd8579d55fc95d9af2a94f892b7c2dc2d88e9`）から、200四辺形の角231点と拘束だけを抜き出したデータである。Apache-2.0の幾何データであり、661節点の2D変位は含まない。元のQ8の4辺を保つ双一次写像を面内で1/2/4分割し、同じ奥行き1層へTet10で押し出す。

準備した1,200／4,800／19,200要素、2,583／9,963／39,123節点は、正のJacobianと体積220 m³の本体 `validate` に合格した。奥行き分割は全段で同じである。F=1の元材料（E=100000 kPa、ν=0.3、γ=20、c=10 kPa、φ=20°、ψ=0°）を使い、過去の約1.35低減済み比較とは別ケースにする。過去の材料や参照値を書き換えたものではない。

次の求解はユーザーPCで行う。この環境では入力準備と検査だけを実施した。PowerShellの例：

```powershell
python tools/prepare_slope_h_series.py runs/griffiths-inplane
foreach ($level in @('xy1-ny1', 'xy2-ny1', 'xy4-ny1')) {
    $p = "runs/griffiths-inplane/$level"
    target/release/soil-fem-3d.exe validate "$p/mesh.json"
    if ($LASTEXITCODE -ne 0) { throw 'mesh validation failed' }
    target/release/soil-fem-3d.exe solve-stages "$p/mesh.json" "$p/ground.json" "$p/case.json" "$p/result.json"
    if ($LASTEXITCODE -ne 0) { throw 'slope solve incomplete; no convergence or Fs verdict' }
}
python tools/assess_slope_h_series.py runs/griffiths-inplane runs/griffiths-inplane-assessment.json --displacement-overlap
```

評価器は準備時の入力SHA、同じ材料/工程・ソルバー、完成と元の物理ゲート、原200四辺形ごとの体積を検査する。共通天端位置の変位には参照なしRichardsonを適用し、応力は原200四辺形プリズムごとの体積平均として比較する。これは積分点応力そのものの全場L2ではない。`--displacement-overlap` は変位の全交差領域積分を実行するため、大きい系列では時間・メモリーを要する。元2Dデータとの5％判定は行わない。奥行き2層の別系列は新しい出力先で `--depth-divisions 2` を指定し、1層の系列へ混ぜない。

## 記録と残件

実測JSONは [2026-11-07の記録](../interim/2026-11-07-plastic-spatial-convergence/README.ja.md) に保存する。Rust本体、従来の検証コード、依存定義は変更していない。全既存試験と新規試験の件数、Windows CIの結果は同記録に追記する。

未達は、細分化で次数が維持されるか、奥行き差、非一意性・局所化、押出しの応力と振動する掘削モニター、Griffiths斜面の実h系列求解・元2D変位5％、正確なFs、既存のMCC掘削274〜14,909要素の空間差、構造・水・JOINTとの統合、一般3D空間収束、P6全体である。

## 奥行き差と追加診断の追補

元の数値と判定は保持し、[奥行き別系列・押出しの追加診断・MCC記録の分離](DEPTH_SPATIAL_DIAGNOSTICS.ja.md)を別に追加した。固定奥行き2層の押出し面内系列も新しい出力先に保存し、元の1層系列へ混ぜない。一般3D空間収束とP6は引き続き未達である。
