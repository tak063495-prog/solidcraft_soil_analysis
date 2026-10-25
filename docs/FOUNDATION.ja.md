# 多層地盤・K0・MCC基盤の操作と仕様

2026-10-23実装。SolidCraftの名前付きボディから、整合する多層体積メッシュと地圧初期状態を作る。本書は初期化と材料点APIの仕様。続いてMCCの全体増分釣合い・Newton/切戻しも実装したため、施工計算は[施工解析仕様](STAGED_FEM.ja.md)を参照する。

## 二層地盤を実行する

```powershell
cd D:\Work\soil-fem-3d
. .\scripts\Enter-Dev.ps1
.\scripts\Invoke-Ground.ps1
```

`examples/layers.json`はSolidCraftのコマンドで、LowerClay（Z=0–5m）とUpperClay（Z=5–10m）を別ボディとして作る。水平範囲はともにX,Y=0–10m。`examples/ground.json`がボディ名と材料ID、単位体積重量、K0、圧密状態、MCC物性を指定する。

保存済みモデルや起動中モデルでも同じ経路を使用できる。

```powershell
.\scripts\Invoke-Ground.ps1 -Design 'D:\path\layers.solidcraft' -GroundFile 'D:\path\ground.json'
# 二層サンプルをアプリで開く（7879が空いている場合）
.\scripts\Start-SolidCraft.ps1 -Design tools\bootstrap-layers.solidcraft -Port 7879
# 起動中モデルを127.0.0.1:7879で取得
.\scripts\Invoke-Ground.ps1 -Live -Connect 127.0.0.1:7879
# 細かいメッシュ
.\scripts\Invoke-Ground.ps1 -MeshSizeM 1.25
```

ボディを追加・削除・改名したら対応表も更新する。未割当ボディを黙って無視しない。SolidCraftの表示材料やSTEPのボディ順から解析材料を推定しない。解析中のCAD変更は取得前後のinspect比較で検出し、途中で変更されたスナップショットを採用しない。

各ボディを別STEPに書き出し、体積とSHA-256を保存する。Gmshのfragment祖先対応から材料IDを追跡し、接触面を節点共有で一括メッシュ化する。正体積の重複、STEPのハッシュ違い、消失ボディ、メッシュ/ボディ体積違いを拒否する。接触するだけの層境界は重複体積として扱わない。

## 初期状態の定義

現在の`horizontal_k0`は、同じ矩形水平範囲を持ち、隙間・重複なく積み重なる水平直方体層だけを受け付ける。一般斜面へK0を代入する機能ではない。各Tet4の重心で全鉛直応力を層別全重量から積算し、既知の静水圧を差し引いて有効鉛直応力を求める。水平2方向の有効応力はそれぞれK0x、K0y倍。水中重量への二重置換はしない。

側面は法線方向ローラー、底面は鉛直拘束。上面の指定上載圧は三角形面から整合節点荷重へ配分する。処方応力から全体内力を組み立て、自由自由度の真の力残差と支持反力を計算する。不釣合いの初期応力をゼロ変位として受理せず、入力エラー/不釣合いとして停止する。この経路では応力を緩和する全体反復は行わない。

水位は未指定で乾燥、または上面/水平層境界/領域外に指定する。重心1点の応力積分を使うため、層内部を横切る水位は現段階では拒否する。地表より高い外水位も外水圧積分を未実装のため拒否する。計算する水圧は指定静水圧で、浸透解ではない。

### OCR・pc・間隙比

`consolidation.definition`で定義を必ず指定する。

|指定|意味|
|---|---|
|`mcc_pc_over_mean_effective` + `ratio`|MCCの等方pc/現在平均有効圧力。鉛直OCRとは別の比|
|`vertical_ocr` + `ratio` + `preconsolidation_pressure_kpa`|鉛直OCRと、別途同定したMCCのpc。鉛直OCRだけではpcを決めず入力を拒否|
|`explicit_pc` + `preconsolidation_pressure_kpa`|MCCのpcを直接指定。鉛直OCRを推定しない|

例の比2は**MCCのpc/p**である。K0による偏差応力があるので、この比を1にすると降伏面の外になる場合がある。その場合pcを自動的に引き上げず拒否する。

間隙比は`void_ratio.mode=from_mcc_reference`で次式から積分点別に設定する。

`1+e = vλ − λ ln(pc/pref) + κ ln(pc/p)`

λ・κはe（またはv）対自然対数の傾き、prefとpcとpはkPa。入力がCc/Csなら利用者が対数底を確認して変換する。`explicit` + `value`も指定できるが、同じ基準線/pc/pに整合しなければ拒否する。深さでpが変わる層に、一定のeと一定のpc/pを独立に与えられるとは仮定しない。出力にe、間隙率e/(1+e)、pc/p、鉛直OCR（指定時のみ）を区別して保存する。

## MCC材料点

有効応力とひずみは**圧縮正**。成分順はxx,yy,zz,xy,yz,zx、ひずみのせん断成分は工学せん断ひずみ。従来弾性ソルバーの引張正とAPIを区別する。p=tr(σ′)/3、q=√(3/2 s:s)、降伏関数は`q² + M² p(p−pc)`。関連流れ則、K=vp/κ、指定νによるG、塑性体積変形に伴うpc発展を用いる。[公開モデル資料](https://docs.itascacg.com/itasca920/common/models/camclay/doc/modelcamclay.html)を参考に独立実装した。

サブステップ開始時のK・Gと硬化係数v/(λ−κ)を固定し、p,q,ln(pc),塑性乗数の4未知量を後退Eulerで戻す。局所Newtonとラインサーチを使い、不成功時は確定状態からサブステップ数を倍増して再試行する。vは各体積ひずみ増分で指数更新する。有限増分では近似誤差があるため、最大サブステップの縮小と別の速度式積分で精度を照合する。これは全体大変形FEMの実装を意味しない。

M>0、λ>κ>0、−1<ν<0.5、p>0、pc>0、e>0、降伏面の許容値を検査する。圧力のゼロ切りや、降伏面を無断で拡大する処理はない。状態には有効応力・e・pc・6成分塑性ひずみを保持する。`integrate`は元状態を変更せず新しい試行状態を返す。JSON読戻し後の次増分が同じ結果になることを試験したが、全体施工の保存再開はまだ実装していない。

```powershell
cargo run --locked --release -- mcc-path examples\mcc-isotropic.json runs\mcc-isotropic-new-result.json
```

出力ファイルは新規のみ。上記を繰り返す場合は新しいファイル名を指定する。

## 出力と検証

`Invoke-Ground.ps1`は、ボディ別STEP/対応表、mesh.json/MSH、初期状態JSON/VTU、独立検算JSONを新しいrunへ保存する。VTUはセル値の材料ID・領域ID・全/有効応力・水圧・e・間隙率・pc等を含む。検算成功時に`accepted-run.json`へ各ファイルのハッシュを記録する。CAD/材料を変更した結果の自動再開機能ではない。

```powershell
# Rust試験、CAD連携、初期状態検算、独立MCC参照、CAD異常系を一括確認
.\scripts\Test-Foundation.ps1
```

最新の数値と検証範囲は[実装検証記録](FOUNDATION_VERIFICATION.ja.md)を参照。

## 今回の実装境界

この基盤に加え、現在は[小モデルのMCC施工解析](STAGED_FEM.ja.md)を実装済み。Tet10、一般斜面の地圧平衡、局所掘削、時間圧密、SRM・構造部材・履歴JOINT・浸透は後続工程である。本書の非排水試験は材料点の体積一定ひずみ経路であり、連成水圧を解いた全体非排水解析ではない。
