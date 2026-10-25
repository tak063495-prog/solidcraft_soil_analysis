# 環境構築と接続検証の記録

2026-10-23、日本時間。作業先：`D:\Work\soil-fem-3d`。

## 完了した項目

- Rust 1.99.1 / cargo 1.99.1 / rustfmt / clippyを設置し、プロジェクトのtoolchainファイルで固定。既存のデフォルトRust 1.85.0は維持。
- 既存のVisual Studio 2022 Build ToolsでWindows x64バイナリをビルド。
- SolidCraft公式v0.20.0のportable ZIPを取得し、公式リリースのSHA-256と照合。
  `3375af3534a3a252bafc12b71d85da86d86e6293922b3dc86975aa45843734a8`
- 専用Python 3.13.14 venvへGmsh 4.15.2、meshio 5.3.5、numpy 2.5.4を設置。`pip check`成功。
- SolidCraftと2D版のソースをsources.lock.jsonのコミットへ固定。両参照checkoutの変更なしを確認。
- 新規プロジェクトをGit初期化。生成物・venv・外部checkout・配布バイナリを追跡対象外に設定。リモートへの公開は行っていない。
- `Setup.ps1 -WithReferenceSources`を既存環境で再実行し、正常終了。空の別PCに対する実機検証は未実施。

## Rustの数値確認

`cargo fmt --check`、`cargo clippy --locked --all-targets -- -D warnings`、`cargo test --locked`成功。3試験：

1. Tet4の一定ひずみ、剛性対称性、並進3・回転3の剛体モード。
2. 反転要素、不正E、ν=0.5/−1等の入力拒否。
3. 底面3節点固定の単一四面体について、解析解の鉛直変位−1/240m、支持反力20/6kNと照合。

これらは接続用弾性実装の試験であり、弾塑性材料則の検証ではない。

## CADから解析までの通し実行

|項目|10m角ブロック・起動中SolidCraftから取得|斜面サンプル・保存モデル経由|
|---|---:|---:|
|run|`runs/20261023-195703-a020f354`|`runs/20261023-200211-10f0a414`|
|経路|MCP→起動中アプリ→STEP|SolidCraft JSON→.solidcraft/STEP|
|寸法範囲 [m]|0–10, 0–10, 0–10|X:0–30, Y:0–20, Z:0–20|
|節点数|140|300|
|Tet4数|373|936|
|体積 [m³]|1000|9000|
|全自重 [kN]|20000|180000|
|CG反復|94|122|
|元行列の真相対残差|8.476×10⁻¹⁰|9.363×10⁻¹⁰|
|反力＋自重の合力誤差/全重量|2.943×10⁻¹¹|3.192×10⁻¹¹|
|最大変位 [m]|0.0959184|0.435922|
|VTU読戻し・変位/要素一致|成功|成功|

E=10000kPa、ν=0.3、γ=20kN/m³。底面のXYZを固定した自重試験。地盤定数や境界条件はサンプル値であり、これらの変位を実地盤の予測として使うものではない。斜面は一般3Dソルバーへの接続形状で、奥行き方向の断面形状は一定。

`tools/verify_run.py`はRustとは別のnumpy行列式で体積を再計算し、既知の寸法範囲、重量、反力、VTUの読戻しを確認する。

```powershell
.\.venv\Scripts\python.exe tools\verify_run.py runs\20261023-195703-a020f354 --volume-m3 1000 --expected-bounds-m 0 0 0 10 10 10
.\.venv\Scripts\python.exe tools\verify_run.py runs\20261023-200211-10f0a414 --volume-m3 9000 --expected-bounds-m 0 0 0 30 20 20
```

各runの`verification.json`に数値を保存。斜面の`cad-slope.png`はSolidCraftのGPU描画で取得・目視確認済み。VTUの配列は読戻して確認したが、ParaViewによる表示試験は未実施。

## 作業中に修正した点

初期の斜面サンプルはSolidCraftのXZ平面のローカル第2軸を誤ってZ正方向と解釈し、上下が反転していた。体積は一致していても支持領域が異なるため、この試行は採用しない。スケッチ座標を修正し、Z=0–20m、底面Z=0を独立確認して再計算した。旧run `20261023-200046-6425f1e4`は`QA_REJECTED.txt`で識別している。

## 起動中の連携

検証時に`Start-SolidCraft.ps1`でアプリを起動し、127.0.0.1:7878の認証付きローカル制御へ`solidcraft-cli mcp --connect`で接続した。アプリはモデル作成に使えるよう起動した状態にしている。閉じた場合は同スクリプトで再起動できる。Codex全体のMCP設定は変更しておらず、プロジェクトのブリッジ経由で接続する。

## 未実装・今回の検証対象外

3D MC、非線形釣合い、SRM、施工、水、JOINT、多ボディ/材料割当、面条件のGUI、SolidCraft内の結果コンターは未実装。外部モデル一般の交差・重複体積の完全検査、巨大メッシュ、特異系すべての診断、配布パッケージのライセンス監査も今回の接続確認には含めない。これらの実装順と受入条件は[計画書](PLAN.ja.md)に記載した。
