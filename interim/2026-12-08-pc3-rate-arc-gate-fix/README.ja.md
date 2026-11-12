# PC-3 rate-arc ゲートの修正と PC-4

PR #13、`codex/plastic-spatial-convergence`、起点 `6b746fa05eafaef626d048270355fb512ab8c522`。フォルダ名の2026-12-08は依頼ラベルで、作業日は2026-11-13（JST）。**NC-01・P6は未達。`passed=false`、`nc01_complete=false`を維持し、PRはマージしない。**

## PC-3の停止記録（ユーザーの転記）

原PCファイルはクラウドから取得していない。以下はチャットで提供されたPC出力の転記であり、今回のクラウド新規計算・原出力監査とは区別する。

- worktree：`D:\Work\scsa_nc_arc_PC3`、HEAD `8342762bd714e33523ce5fef772572699fc21881`、出力 `runs\NEW-nc-fold-confirmation-PC-3`。
- 2026-11-12 20:52 JST、公式 `interim/2026-12-05-nc-independent-arc-residual/rerun-pc.ps1` のrate-arc400直後の旧ゲートで停止。`pc3-failure.json` は `passed=false`、`nc01_complete=false`、`stop_reason="rate-arc400 不足/残差超過/停止。失敗を保存"`。
- `rate-arc400/report.json` は `equilibrium_samples=128`、`stop_reason=null`、`maximum_bordered_true_linear_residual=5.4291880423103136e-11`、`passed=false`。`samples.json` は128件。
- runtime / build / focused-tests / reference400 / rate-arc400 はすべてrc=0。independent400-coarse/fine/material8と監査、800系、§4.1、§5は未実行。

## 原因と仮説の切り分け

旧コードは `Read-Json` が `Get-Content -Raw -Encoding utf8 | ConvertFrom-Json` をそのまま出力し、呼出し側が `$Samples = @(Read-Json ...); $Samples.Count` としていた。**JSON配列を単一の配列オブジェクトとして出力するシェルでは、外側の `@(...)` の件数は1である。内側の実際の128件を数えず、`$Samples.Count -ne 128` だけで停止する。** JSONの中身が二重配列なのとは異なる。

PowerShell 7.4.6の実測では列挙され、旧式でも128件となった。このため、提示された値だけでPowerShell 7上の停止を再現したとはしない。Windows PowerShell 5.1と7の実機比較には、[同じ旧式と新ゲートを検査するスクリプト](test-rate-arc-gate.ps1)と[専用CI](../../.github/workflows/pc-rate-arc-gate.yml)を追加した。PC-3自身のPowerShell版・生JSONを未取得なので、PC-3の実環境との最終照合はPC側の結果で行う。Windows CIの実測結果は下の検証記録へ追記する。

|候補|コード・試験から確認したこと|
|---|---|
|配列列挙と件数|旧式はシェルの出力列挙に依存する。7.4.6では128。新ゲートはPythonでトップレベル配列を直接数える|
|生成JSONの形式|`independent_arc.py` は `frames=[]` に辞書をappendし `json.dumps` で保存。トップレベル配列と辞書行。実生成1点の試験でも確認|
|stop_reason|停止なしはPython `None` → JSON null。文字列 `"null"`、空文字、停止文字列は拒否。キー欠落も拒否|
|passed|生成器はPython `False` → JSON boolean false。欠落・null・数値0・文字列 `"false"` は受入にしない|
|数値比較|件数はJSON整数128、残差は有限のJSON数値で `0 <= r <= 1e-10`。提示値は合格。数値文字列・bool・非有限値は拒否|
|Read-Json|JSONテキストを書き換えずPowerShellオブジェクトとして出力する実装。配列列挙の依存を新しいrateゲートから除いた|

実生成試験は元のSHA固定 `reference400-fixture.zip` から既存生成器を**1点だけ**実行する。128点の数値計算は行わない。再利用検査の試験に使う128行はその1点を複製した**合成fixture**であり、PC-3や独立経路の数値的証拠ではない。

## 修正と限定受入

[`tools/nc_pc_gate.py`](../../tools/nc_pc_gate.py)を両シェル共通の検査器とし、公式[rerun-pc.ps1](../2026-12-05-nc-independent-arc-residual/rerun-pc.ps1)の400/800ゲートから呼び出す。許容値・参照値・数値求解ソースは変更しない。

- 合格条件は128点、停止理由が明示的null、元の拡張系残差が1e-10以下、booleanの `passed=false`。空配列・欠落キー・件数不足・入れ子配列・不正な型を拒否する。`accepted=true` はこの限定ゲートの通過だけを表し、診断出力の `passed/nc01_complete/p6_complete` はすべてfalse。
- `rate-arc*-gate.json` に元のreport値、samplesの型・件数、全不合格条件 `failures` を保存する。外側の `pc3-failure.json` / `pc4-failure.json` に実行段階、診断パス・内容、`failed_conditions` を含める。ネイティブ呼出しのrc・引数・時刻も従来どおり残す。
- 新規出力先だけを許可し、検査JSON・失敗JSONは排他的新規作成。既存失敗を消さず、再実行で上書きしない。PowerShell版も新出力へ保存する。
- `-ResumeFrom` はPC-3のreference400/rate-arc400だけを読み取り参照する。元HEAD・rc・引数、400/RK4前史のexcavate_upper起点、元入力、128サンプルの形状・有限値・制御、旧HEADとの計算ソース/入力同一性を検査する。旧HEADのtools Python全件・対象arc Python・examples JSONを保守的に照合するため、無関係な変更でも停止する場合がある。検証失敗は再利用せず停止する。
- 参照した400系全ファイル、元HEAD・コマンド・PC-3失敗記録のSHA256をPC-4へ保存し、正常終了前に再照合する。大きな線形系ファイルはストリームでhash化する。元ディレクトリ内への出力は拒否する。

これらのソフトウェア検査と小規模fixture試験を限定受入とする。PC-4の原出力検査・残りの実行手順は[REPRODUCE.ja.md](REPRODUCE.ja.md)に記載した。

## 未達を維持する項目

**PC-4は未実行。** independent400の3本と監査、reference800/rate-arc800/independent800の3本と監査、§4.1の直接64対32+32、§5のnative対独立比較はPCで実行する。全場・弱モード・荷重・刻み/前史/材料内4対8の感度、peak一致基準の未定義、全上部撤去、全解放後凍結、本番完了までの再開は残件。ゲート通過から局所fold確定・物理安定限界・NC-01/P6完了へ昇格させない。

案2・正則化・粘性材料・p'下限・成分別解放率・追加拘束・不釣合い力再配分・疑似過渡/動的緩和は実装しない。過去の失敗・未実行・不合格記録は保持する。

## 検証記録

最終Linux検査は **41 passed / 0 failed / 0 skipped**（新規21件＋既存20件）。[validation.json](validation.json) と `validation/` にコマンド、rc、各件数、原ログ、実行前後のソースSHAを保存した。Rust数値ソースは変更せず、大規模計算・全件Rust/Python試験は今回実行していない。

|対象|passed|failed|skipped|
|---|---:|---:|---:|
|新規 `test_nc_pc_gate.py`|15|0|0|
|新規 `test_nc_pc_gate_schema.py`|6|0|0|
|既存 `test_nc_independent_arc_residual.py`|6|0|0|
|既存 `test_nc_fold_stability.py`|14|0|0|

PowerShell 7.4.6の[実測JSON](powershell-7-result.json)は10ケース成功・失敗0・skip0、旧式128件/直接代入128件。シェルなしの別実行も **14 passed / 0 failed / 1 skipped** を確認した。これらの再実行や実スクリプト試験内で呼び直す回帰を、41件の一意試験数へ加算しない。

初回の既存fold検査は **13 passed / 1 error / 0 skipped**。疎配置で元 `native-evidence.zip` が未配置だったためで、固定HEADの元blobを復元後に14件全件を再実行して成功した。`validation/initial-fold-fixture-error.log` を保持し、初回が成功だったとは記録しない。

Windows PowerShell 5.1/7の専用CIは初回push時点では未確認。実機結果を確認後に追記し、旧先端のCI成功を今回の証拠へ転用しない。
