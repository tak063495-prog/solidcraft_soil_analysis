# Windows CI の旧ワークフロー SHA 検査

PR #13 の先端 `20a69407bee17bdeebebdb7e10025f9e91a9f9f3` に対する Windows の python-core-0 は 475 件中 473 成功・2 失敗。A-03 と B-03 の保存試験がともに同じ SHA256 の差を報告した。[run 44316604987 のログ](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/44316604987/job/134844266293)

`18d54b1` が core-0/core-1 用の PC 診断アーカイブ選択行を追加した。親 `520cecf` の SHA256 は `cc0beddc1769051a2f68d33dee3057ec19c07df77bad8de5a982ea637346e2d5`、追加後は `2207dbc1e9ca37e7c30ef2f41f327faebc427498f96d6a4407ada4fc59620287`。その後 `20a6940` まで同ファイルの差分はない。以前の Linux 集計は両試験に初回実行の成功を採用し、再開後に保存したソース指紋は追加後の SHA だった。今回、先端の同一ファイルで両試験の失敗を再現した。

旧マニフェストを使い続けるため、追加行の字下げ、引用符、条件、パス、末尾 LF、出現回数 1 と次ステップ直前の位置を検査する。照合した 1 行だけを除いた全バイトを旧 SHA256 と比較する。他のファイルは原本全体を照合する。`tests`、`tools`、Rust 試験、過去の保存マニフェストとその読込箇所を検索した。現在のワークフローをこの旧 SHA で拘束する通常試験は A-03・B-03 の 2 件。過去のソース指紋やアーカイブ再生は当時の記録を参照する。

P6-01 の追加出典はこのフォルダの acceptance.json。固定コミット、Git blob ID、SHA256 と各項目の型付き false を既存監査器で照合する。

## 全件試験の実測

対象は `006d0d114477e6913c2a8eb845294a8350692df4`。Rust は全バイナリの一覧と実行 ID を、Python は 7 群の割当と実行 ID を照合した。

| Linux | 対象 | 成功 | 失敗・error | skip | アーカイブ不足で未実行 |
|---|---:|---:|---:|---:|---:|
| Rust 全ターゲット | 467 | 467 | 0 | 0 | 0 |
| Rust doc-test | 0 | 0 | 0 | 0 | 0 |
| Python（tests・tools） | 981 | 979 | 0 | 2 | 0 |

Python の skip は PowerShell を呼ぶ 2 件で、core-0 と core-1 に各 1 件。両群のプロセス終了値は既存ランナーの skip 判定により 1。各 ID と実出力は [validation.json](validation.json) と [Linux の実行ログ](linux-results.zip) に保存した。ZIP には取得した 1,076 パスの blob ID とバイト数も収録した。

同じ PR 先端に対する Windows FEM は成功。core-0 は 475 件、core-1 は 482 件すべて成功し、旧 SHA 検査 2 件と PowerShell ゲート 2 件の個別成功を照合した。[Windows FEM 44396740515](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/44396740515)・[ゲートと旧試験の実ログ抜粋](windows-core.log)・[Windows PowerShell 互換性 44396740610](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/44396740610)

## 保存アーカイブの取得

PR の履歴は shallow ではなく、チェックアウトは `blob:none` の partial clone。失敗した `007e4c9b…` は先端 tree にある 14,361,601 バイトの通常の blob で、対象パスの LFS filter は未指定だった。コミット単独の fetch は成功したが、通常の blob fetch は filter の有無どちらでも終了 1。同じ OID を `--refetch --no-filter` で再要求すると取得できた。別の失敗 OID では個別の通常 fetch も成功したため、全件に共通するサーバー内部の原因までは確定していない。

GitHub の固定コミットの raw URL から試験用に選んだ 1,076 ファイルを取得し、すべて先端 tree の blob ID と照合した。前回記録された 650 の依存パスをすべて含む。取得時の比較は [archive-diagnosis.log](archive-diagnosis.log)、修正前の Linux 再現は [linux-before.log](linux-before.log) に保存した。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
