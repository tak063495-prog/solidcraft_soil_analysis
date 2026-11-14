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
P6-03 の検査もそのままで、P6-01 への連鎖が解消したことを Linux と Windows の実テストで確認した。

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
|通常 Python 全件|**676 固有 ID すべて成功**。core-0: 321 / core-1: 347 / mean: 4 / spatial: 4。最終 failure/error/skip/expected failure/unexpected success はすべて 0|
|旧先端 5e3d191 の Windows|run 43325715202 は failure。両 core 群が failure、他 6 群は success、集約は failure。PC compatibility 43325715212 は success|
|修正後 Windows Python|core-0 321 / core-1 347 / mean 4 / spatial 4、計 676 件を実行し、4 群すべて success。両固定監査の `... ok` をログでも照合|
|今回の修正の Windows|**08b58fef4727659ad91e2fafaa087aef53e7ee52** の [FEM 43330893373](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43330893373) は完了・success。8 検証群と集約の全 9 job が success。[PC compatibility 43330893601](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/43330893601) も両シェル success|

Windows は上記の **同一先端 SHA** と run 結論を照合し、FEM 全 9 job と compatibility 全 2 job の完了・success を確認した。
[Windows job 一覧](validation/windows-08b58fe-jobs.tsv) と [Python 件数・該当監査のログ抜粋](validation/windows-08b58fe-python-excerpts.log) を保存した。
キャンセル・時間切れ・実行中を成功へ数えていない。旧 5e3d191 の failure は [旧 job 一覧](validation/windows-5e3d191-jobs.tsv) にも保持する。

この追補は検証済みコード 08b58fe に結果文書・ログを追加するもので、数値・監査・試験ソースと台帳登録は変更しない。
結果記録コミットによって発火する新しい CI は、[PR #13 Checks](https://github.com/tak063495-prog/solidcraft_soil_analysis/pull/13/checks) で別に確認し、
作業完了時の返信でも先端と結果を明示する。この文書内の成功記録の対象 SHA を後続の文書コミットへ付け替えない。

全件試験前に固定した **3,919 ファイルは終了後もすべて同一**。
最終 4 群は全て終了 0。各群の件数・終了コード・実測時間・ログ SHA は [linux-summary.tsv](validation/linux-summary.tsv)、
重複のない全 ID は [linux-test-ids.tsv](validation/linux-test-ids.tsv)、生ログは `validation/linux-python-*.log` に保存した。
該当 57 + 19 件は通常集合にも含まれるため、676 件へ加算しない。
初回 core-0 では疎配置から漏れた旧 cut-volume `evidence.zip` による準備エラーを検出し、
元コミットの blob と完全一致するファイルを復元した。初回 core-0 の 321 件は 320 成功・準備エラー 1・skip 0、終了 1。
同群全 321 件を元の CI runner で再実行し、全成功・終了 0 を確認した。
初回の全群 676 件と同群再実行 321 件の method 試行は **997 件**。
最終結果は 676 固有 ID に対して判定し、再実行や該当試験を水増ししない。
初回ログ・初回集約の終了 1 は `linux-core-0-first-attempt.log` / `linux-initial-processes.log` に保持した。
これはテスト期待値や証拠の変更ではない。

Git tree では起点の旧 6,933 パスをすべて保持し、監査器・試験以外の 6,931 blob ID と実行属性は不変、削除は 0 件。
起点 6,933 追跡ファイルのうち、配置した旧 5,898 ファイルはバイト・実行属性が一致した。
既存変更は監査器と試験の 2 ファイルのみ。旧 JSON 2,687 件はすべて同一バイトで、
省略したアーカイブ等 1,033 パスの中身は今回検証したとはしない。
内訳は `validation/preservation.tsv` と `validation/sparse-omissions.tsv` に保存する。

## 保持した限界

数値ソルバー、許容値、参照値、既存の数値証拠・過去の失敗は変更していない。
NC-01 の全撤去、B-02 の 768 Tet10 の元残差基準、一般空間収束、P3b/P4、P6 全体を
CI の成功から合格へ昇格させない。キャンセル、時間切れ、実行中は成功に数えない。
