# 混合MCCと梁・板の連成：限定した独立検証

2026-10-26。Tet10変位・材料境界で分離したP1平均有効応力のMCCへ、線形3D梁とCST/DKT板を接続した。独立GP体積形式とP1投影体積形式を区別して検証する。これは固定した二層試験体の受入であり、実CADの塑性山留め、空間収束、低圧域、水圧/JOINT、native部材表示やP6全体の受入ではない。

## 変更

- 混合方程式の変位ブロックへ部材の内力と接線を加える。疎行列には、部材の縮約で生じる離れた土節点間の結合も登録する。
- 設置・撤去の解放荷重と確定した部材変位・回転・設置基準を、切戻し・チェックポイントへ保存する。
- 結果形式を `soil-fem-mixed-mcc-stage-result/2`、`soil-fem-projected-mcc-stage-result/2` とする。部材のない従来の `/1` は継続し、投影型の材料応力と釣合いに使う仕事共役応力、材料体積と運動学的体積を区別する。
- VTU/CSVと解析パッケージ監査は新形式に対応する。nativeビューアの既存形式へ読み替えることはしない。

## 試験体と独立解

`tools/prepare_mcc_structures.py` は既存の二層・12 Tet10・324積分点の柱を使う。K0、平均有効応力を基準にしたOCR、MCC基準線の間隙比から初期状態を作る。140 kPaプレロード、変形後の梁・板設置、部材載荷、180 kPa硬化、100 kPa除荷、上層0.5 m³掘削、部材除荷、撤去の8工程・9状態で比較する。梁と片持ち板の土への接続は残存する中間面に置く。掘削後も100 kPaの圧力を新しい地表に加える試験であり、完全無載荷の自由掘削底面への一般化はしない。

`mcc_structural_reference.py` と `run_mcc_structural_reference.py` は、元入力だけから初期状態・土の塑性履歴・設置基準を作る。MCCの解析的な速度形接線と事象検出付き中点積分を使い、梁・板の自由並進、回転、端部解放を混合KKT系に残す。Rustの縮約行列、解いた変位、塑性履歴、解放ベクトルは参照解の生成に使わない。梁は断面柔性と仮想仕事、板は独立した多項式DKT拘束で組み立てる。

参照解は100/200分割、Rustは施工刻み1/8→1/32、材料最大刻み1e−4→2.5e−5で比較する。比較ツールは、参照解の全場差1e−4以下、両Rust解の全場・部材相対差1％以下、Rust細分化で主要5場の最大差が半分以下になることを要求する。圧力下限、応力の投影修正、人工剛性、残差基準の緩和は導入しない。

## 完了した照合

|方式|最大変位相対差：粗→細|最大応力相対差：粗→細|最大塑性ひずみ相対差：粗→細|
|---|---:|---:|---:|
|P1投影体積|0.35015％→0.08709％|0.03374％→0.00845％|0.27376％→0.06605％|
|独立GP体積|0.34983％→0.08737％|0.03251％→0.00811％|0.27491％→0.06503％|

pc、間隙比、全節点変位、全GP応力/塑性履歴、部材並進/回転、設置基準、梁端力、板断面力/モーメント、設置/撤去/掘削の解放力が各ゲートを通過した。硬化後は324積分点すべてに塑性履歴を持ち、上層撤去後の残存側は162点。参照解の力・体積残差も1e−8以内に保つ。

今回のRust版で全249件、Python全94件とfmt/clippyが合格した。固有の接線検査2件、回転した3D梁/板/端部解放を含む独立全体系の検査2件、構造連成出力の正常/改変拒否/重力平衡基底荷重の検査6件、元ソース保存の検査3件を含む。重力平衡の基底荷重を重複計上しない回帰検査も追加した。

両方式で各96受理増分・4意図棄却、細側は各384受理増分・4意図棄却。粗側の各8工程、計16地点で途中保存→別プロセス再開し、全状態・フレーム・部材履歴・監査・増分経路が完全一致した。4計算36フレームのVTUと生値CSVも監査した。部材モーメント、設置基準、解放力、結果形式、必要な報告を改変した拒否試験も通過する。

記録：

- [P1体積の全場・部材比較](../runs/mcc-structures-projected-20261026-v1/comparison.json)
- [GP体積の全場・部材比較](../runs/mcc-structures-gp-20261026-v1/comparison.json)
- [P1体積8地点再開](../runs/mcc-structures-projected-20261026-v1/checkpoints/verification.json)
- [GP体積8地点再開](../runs/mcc-structures-gp-20261026-v1/checkpoints/verification.json)
- 固定実行ファイル：`runs/mcc-structures-20261026-v1/solver/soil-fem-3d.exe`、SHA-256 `7738b0e38f18a19c5872ed663e7d7d0264f5759bd8947c662ccb48f41f66cb2f`。

## 再現

```powershell
cargo build --locked --release
python tools/prepare_mcc_structures.py runs/new-mcc-structures
python tools/run_mcc_structural_reference.py runs/new-mcc-structures/ground/mesh.json runs/new-mcc-structures/ground/ground.json runs/new-mcc-structures/stages.json runs/new-mcc-structures/reference100.json --steps 100
target/release/soil-fem-3d.exe solve-stages runs/new-mcc-structures/ground/mesh.json runs/new-mcc-structures/ground/ground.json runs/new-mcc-structures/stages.json runs/new-mcc-structures/stage-results.json
python tools/export_stages.py runs/new-mcc-structures/ground/mesh.json runs/new-mcc-structures/stage-results.json runs/new-mcc-structures/export
python tools/verify_stage_package.py runs/new-mcc-structures/ground/mesh.json runs/new-mcc-structures/ground/ground.json runs/new-mcc-structures/stage-results.json runs/new-mcc-structures/export
```

GP体積には準備時に `--gp-volume` を付ける。比較用の細側は同じ物理入力で数値刻みだけ変更し、参照200分割も別ファイルに保存する。既存の完了記録は上書きしない。

参照結果は実行時の入力・ソースSHA-256を保持する。開発後に元のソースが変わった場合、`verify_mcc_structural_study.py --source-archive` に実行時の `source-snapshot.zip` を指定すると、同一ハッシュの元ツールだけを照合できる。変更されたメッシュ・地盤・工程入力を過去版ZIPで代替することは拒否する。

## 残件

共有面全体を板へ接続する実CADの塑性掘削、土/板の変位跡の整合と空間収束、部材を失う掘削の方針、混合MCCの水/接触/アンカー/土の誕生、nativeの入力と部材結果表示を引き続き検証・実装する。`shared_interface` の実CAD弾性掘削とは別の受入範囲として保持する。
