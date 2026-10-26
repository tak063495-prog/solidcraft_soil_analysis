# MCC中点積分と実CAD壁の途中成果

選択式の二次精度MCCを全体非線形FEMへ接続しました。K0/OCR/間隙比を保持し、弾性経路の厳密積分、塑性中点・初回降伏を含む整合接線、切戻しを実装しています。既定の後退Euler方式は保持します。

固定版のRust263件・Python116件、fmt/clippy/CLIビルド、独立Radauの3経路が合格しました。全分割の状態誤差は128/256/512分割でほぼ4分の1ずつ減少します。実CAD壁の全8工程を独立RK4の100/200分割と比較し、全場・部材・保持力・刻み細分化の原基準が合格しました。細かい解析の全節点変位の最大相対差は0.003743％です。

掘削0.03125・除荷0.0625・壁撤去0.0625の3地点から別プロセス再開し、実測秒数を除く全結果が連続実行と一致しました。撤去1134積分点の履歴と3試行の切戻しを確認し、積分方式を変更した再開は拒否しました。

粗い施工刻み1/8の微小回転の不合格、反復許容値だけを100倍厳しくした対照も残します。施工刻み1/16・1/32で原基準を通過しています。独立速度形計算と本体の材料積分は異なる方法であり、同一メッシュ上の比較を空間収束へ読み替えません。別の局所掘削地盤の参照収束、1792要素の参照、新規PC配布、他のP0〜P6要件は継続中です。完成版ではありません。

## 同じ入力での再実行

リポジトリ先頭から実行します。入力ZIPにはJSONだけを含み、SHA-256確認後、新しい出力先へ展開します。

```powershell
cargo build --locked --release
python interim/2026-10-27-midpoint-mcc/unpack_inputs.py runs/my-midpoint-wall
target/release/soil-fem-3d.exe solve-stages runs/my-midpoint-wall/coarse/ground/mesh.json runs/my-midpoint-wall/coarse/ground/ground.json runs/my-midpoint-wall/coarse/stages.json runs/my-midpoint-wall/coarse/stage-results.json
target/release/soil-fem-3d.exe solve-stages runs/my-midpoint-wall/fine/ground/mesh.json runs/my-midpoint-wall/fine/ground/ground.json runs/my-midpoint-wall/fine/stages.json runs/my-midpoint-wall/fine/stage-results.json
python tools/verify_mcc_midpoint.py --solver target/release/soil-fem-3d.exe --output runs/my-midpoint-material
```

独立壁参照は `tools/run_mcc_structural_reference.py --help` を参照し、同じcoarse入力、`--identity`、`--integration-method rk4` と100/200分割を指定します。`tools/verify_mcc_structural_study.py` で比較します。保存再開は `checkpoint-stages` と `resume-stages` を使い、この入力の受入地点は81・161・177確定増分です。

SolidCraft保守版60ファイルの総合CIと、前回公開bed5e8aのWindows CIも合格しました。新コミットのWindows CIは公開後に別途確認します。検証した実GUI/GPUは従来57ファイルの実行版に対応し、新MCC方式の新規アプリ配布受入を意味しません。コード・正確な入力・試験ログ・検証要約を公開し、実行ファイル、全run生データ、所有トークンは同梱しません。原本の相対パスとSHA-256を各報告へ記録しています。


## 局所掘削地盤の追加比較

初期施工刻み1/256・材料刻み2.5e-5の5工程は、200分割参照に対する全場基準を満たしました。pc最大差は0.0000201 kPa、変位最大差は約0.018 µmです。独立参照100→200分割の変位差約0.0362 µmは上限0.01 µmを超えるため、計算中の400分割で確認するまで総合刻み受入は保留します。`ground-case` に同じメッシュ・地盤・工程入力を保存し、`solve-stages` で再実行できます。壁の入力ZIPとは別のケースです。
