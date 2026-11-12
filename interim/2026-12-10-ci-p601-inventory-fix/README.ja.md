# FEM CI：P6-01 固定ソース台帳の追補

対象は PR #13、`codex/plastic-spatial-convergence`。起点は
`5e3d191b66b04a877728e8b96ae8b027795ba15e`。
保存先の 2026-12-10 は指定されたラベルであり、作業日は 2026-11-13 JST。
**NC-01、B-02 全場収束、P6-01 の数値受入、P6 全体は未達のまま。**
この修正は監査資料の登録であり、数値的な受入を追加しない。

## 原因

P6-01 の `audit_status.py` は `interim` 配下の `acceptance.json`、
`measurements.json`、`validation.json` を再帰検査し、台帳にないソースを
`SOURCE_NOT_INVENTORIED` とする。既存台帳は古い起点 `8342762` の証拠に固定されている。
したがって、新規フォルダの対象 JSON を追加する際には台帳の追補も必要になる。

未登録だったのは次の 3 ファイル。

- `interim/2026-12-08-pc3-rate-arc-gate-fix/validation.json`
- `interim/2026-12-09-b02-tet768-reassembly-residual/measurements.json`
- `interim/2026-12-09-b02-tet768-reassembly-residual/validation.json`

run [43303948088](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43303948088)
の python-core-1 は最初のファイルで失敗、python-core-0 は P6-03 が同じ監査を呼ぶため
`P601_INTEGRITY` で失敗している。ログの失敗箇所と修正前のローカル実行を照合した。
修正前の記録は [validation.json](validation.json) と `validation/before-targeted.log` に保持する。
`validation.json` の `passed=false` は**修正前の記録**であり、最終試験の結果ではない。

## 修正

旧 `ledger.json`、`ledger.ja.md`、引用、行分類、固定起点、過去の失敗記録は変更しない。
新しいソースだけを [inventory.toml](inventory.toml) へ追補する。
上記 3 ファイルは `5e3d191` に固定し、今回の `validation.json` も先行する保存コミット
`24b5d11915f4671f40ec658457a3b00508c01f21` へ固定する。
未コミットの自分自身を参照する循環ハッシュは作らない。

監査器内で追補台帳全体の SHA-256 を固定する。各ソースについて、
存在、SHA-256、Git blob SHA-1、指定した不変コミットの tree 内 blob、
JSON の妥当性、型付き boolean/null の値、根拠のない全体完了を検査する。
旧ソース ID・パスとの衝突、可変 ref、外部 snapshot への付け替え、不正パスは拒否する。
旧ソースには従来どおり `8342762` 固定を適用する。
`SOURCE_NOT_INVENTORIED`、integrity 判定、旧テストの期待値は維持する。
P6-03 の検査もそのままで、P6-01 への連鎖が解消することを確認する。

## 今後、新規 interim フォルダを追加するとき

1. 新規 JSON の内容と `passed`・未達状態をレビューする。既存の証拠や登録を再封印しない。
2. 新規証拠だけを先行コミットへ保存する。この時点ではブランチを push する必要はない。
3. `register_sources.py --commit <40桁SHA> --id <未使用ID> interim/.../validation.json`
   で、そのコミットと作業ツリーのバイト一致を確認し、登録を生成する。
   `measurements.json`、`acceptance.json` と今回のフォルダ自身の JSON も対象に含める。
4. 生成された `[[sources]]` を `inventory.toml` に追記する。`schema` は 1 行だけ保持し、
   既存エントリーのパス・コミット・hash・状態値は書き換えない。
5. 台帳全体の SHA-256 を計算し、P6-01 監査器の `SUPPLEMENTAL_SHA256` を更新する。
   差分をレビューする。登録後に証拠を変更する必要があれば、旧証拠を保持した別パスに保存する。
6. P6-01 / P6-03 の該当試験と全 Python 試験を実行してから、登録とコードをコミットし push する。
   Windows は同じ先端 SHA の全ジョブと集約 job の実際の完了・success を確認する。

具体的なコマンドは [REPRODUCE.ja.md](REPRODUCE.ja.md) を参照。

## 最終検証

|検証|保存時点の実測|
|---|---|
|Linux release build|終了 0|
|P6-01 該当試験|57 成功、failure/error/skip 0|
|P6-03 該当試験|19 成功、failure/error/skip 0|
|通常 Python 全件|676 固有 ID を 4 群へ重複なく割当、実行中。全件成功は未確認|
|旧先端 5e3d191 の Windows|run 43325715202 は failure。両 core 群が failure、他 6 群は success、集約は failure。PC compatibility 43325715212 は success|
|今回の修正の Windows|push 後に同一先端の結果を確認する。現時点では成功未確認|

全件試験前に 3,919 ファイルの指紋を固定した。全件の最終件数・終了コード・生ログは後続の追補で記録する。
初回 core-0 では疎配置から漏れた旧 cut-volume `evidence.zip` による準備エラーを検出し、
元コミットの blob と完全一致するファイルを復元した。初回 321 件は 320 成功・準備エラー 1・skip 0。初回ログを残し、同群 321 件を再実行中。
これはテスト期待値や証拠の変更ではない。

## 保持した限界

数値ソルバー、許容値、参照値、既存の数値証拠・過去の失敗は変更していない。
NC-01 の全撤去、B-02 の 768 Tet10 の元残差基準、一般空間収束、P3b/P4、P6 全体を
CI の成功から合格へ昇格させない。キャンセル、時間切れ、実行中は成功に数えない。
