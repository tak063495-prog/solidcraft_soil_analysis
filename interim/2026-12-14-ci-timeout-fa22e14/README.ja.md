# PR #13: fa22e14 Windows CI の時間切れ対策

対象ブランチは `codex/plastic-spatial-convergence`、起点は
`fa22e14701eaa5014b1d8268e768c1ab513db4a6`。このフォルダの日付は指定ラベル。
PR はマージしない。試験の許容値・参照値・保存値・過去の失敗は変更しない。
NC 関連の数値コード・設定も変更しない。

## 症状と確認した原因

FEM build and verification の run **43446323125**、
job **131956600215** (`python-core-1 on Windows`) は **cancelled**。
集約 `Rust and Python on Windows` は **failure**。成功した CI とは扱わない。

ログ（UTC）では checkout は 14:49:03.262 → 15:26:06.539、約 37 分 3 秒。
最初の `git fetch --filter=blob:none` は約 1 秒で終わり、
`git checkout` の開始 14:49:08.878 から最初の `Updating files` 15:23:55.615
まで約 34 分 47 秒かかる。大きな保存 blob の遅延取得と配置が支配的で、
全履歴の commit/tree 取得だけが 37 分かかったわけではない。
旧ログには Trace2 等がなく、遅延取得のネットワークと runner の展開・ディスク処理の内訳は分離できない。
LFS は false。`Preserve exact archived source bytes` は
checkout 前の `git config --global core.autocrlf false` で、アーカイブ生成・展開の処理ではない。

CLI build は約 3 分 59 秒。Python 開始は 15:30:44.902、取消は 15:33:58.697。
runner の開始 14:48:58.783 から約 44 分 59.9 秒で、定義の
`timeout-minutes: 45` と一致する。取消時のスタックは
`test_c01_path_error_control.GlobalPathErrorSavedEvidence.test_all_history_metrics_are_recomputed_and_forged_acceptance_is_refused`
→ `evidence.py: replay` → `measure` → `study.py: decode` の `KeyboardInterrupt`。
assertion failure ではなく、ジョブの期限到達による中断と判断する。
公開されたログは取消理由を別の専用コードでは示さないため、スタックだけで理由を断定せず、
設定と経過時間を併せた診断である。

同 run の他の 7 検証群は success、checkout は約 5 分 32 秒〜12 分 2 秒。
全ジョブ共通の履歴・アーカイブ容量に加え、転送・ディスク処理の変動がある。
新規 cause archive 359,243,008 bytes だけでは、直前 run からの大幅な遅延全体を説明できない。

比較対象 e8133f0 の run **43427395029**、core-1 job **131892091478** は success。
ジョブ全体が約 22 分、checkout 約 5 分 59 秒、Python は **363 件 / 644.042 秒**。
fa22e14 の追加は cause-diagnosis **8 件**で、全件は 718 → 726 件。
名前順の round-robin のため新しいモジュールの挿入で所属も変わり、
core-1 は 363 → 389 件、core-0 は 347 → 329 件となった。
旧先端の保存 Linux 計測 core-1 379.876 秒、fa22e14 側 625.092 秒の差は約 245 秒だが、
配置変更と実行環境差を含むため、8 件だけの純増とはしない。

## 対策

- C-01 の switch / path / cause のモジュールをそれぞれ専用 Python ジョブへ移す。
  module/class fixture を分割せず、全件を厳密に 1 回ずつ割り当てる。
- 3 つの C-01 archive の合計 **1,917,665,821 bytes** を全ジョブへ無条件に配置するのをやめる。
  switch job は switch の 596,504,815 bytes、path job は path の 961,917,998 bytes、
  cause job は再検証に必要な 3 archive すべてを取得する。JSON/index/source は全群に残す。
  新設する C-01 の 3 job は無関係な zip / tar / gzip / npz を配置せず、
  path replay の retained-failures/prototypes.zip は明示的に戻す。
  追加ジョブによる無関係な大容量 archive の重複ダウンロードを避ける。
  C-01 主 archive の配置量は旧 8 群合計 15,341,326,568 bytes から
  3,476,088,634 bytes へ減る。これは論理ファイル量で、Git 転送圧縮後の実測量ではない。
  path replay は switch の index を読むが、switch archive 自体を必要としない。
- 固定コミット・SHA 検証の台帳と他の保存証拠が参照する履歴を守るため、
  `fetch-depth: 0` と blob 遅延取得を維持する。単純な shallow 化は採用しない。
  既存 NC archive の除外パターン、NC の試験・数値設定は変更しない。
- ジョブ上限は **75 分**にする。37 分の checkout の再発も予算に含め、
  完走に必要な余裕を明示的に確保する。cancelled/timeout は集約で引き続き failure になる。
- 新規 `tools/run_ci_python_evidence_shards.py` が各試験の実測時間も JSON に記録する。
  既存 `tools/run_ci_python_tests.py` は PC-3 の保守的な全 tools ソース照合に含まれるため、
  元のバイトを維持する。NC のゲート・試験・設定を変更して新 runner を除外することはしない。
  replay の中身・判定・保存値は変更しない。

実際の Git sparse checkout を使う回帰試験で、全群の source 保持と必要 archive の配置、
全試験 ID の一意な割当、CRLF で保存した文書の末尾一致を検査する。

## 検証と時間予算

全件の件数・終了状態・実測時間は `validation.json` と `validation/` の生ログに記録する。
37 分 3 秒の checkout、準備約 1 分、CLI build 最大観測約 4 分 32 秒に、
Windows で実測された mean evidence 約 14 分 16 秒を加えても約 57 分で、
75 分まで約 18 分の余裕がある。Rust integration は同 run 最大約 28 分 12 秒で、
同じ checkout 外れ値と準備を足して約 66 分、余裕約 9 分。
これらは予算見積りであり、新先端 Windows の完走保証・成功実測ではない。

## 最終 Linux 検証

Rust **463 件**、Python **730 件**は全て成功。最終の failure/error/未実行/skip/ignore/filtered は 0。
fmt、clippy (`--locked --all-targets -- -D warnings`)、locked release build も終了 0。
Rust の `--list` の 463 件と実行結果を照合した。全 Rust 実行は初回の全ターゲット build を含め
3045.856 秒。69 integration targets、lib/bin、doctest と通常の example build を含む。

| Python 群 | 成功件数 | 実測秒（discovery/fixture を含む） |
|---|---:|---:|
| python-core-0 | 323 | 755.008 |
| python-core-1 | 383 | 671.070 |
| python-mean-evidence | 4 | 1313.089 |
| python-spatial-evidence | 4 | 301.689 |
| python-c01-switch | 4 | 89.092 |
| python-c01-path | 4 | 122.755 |
| python-c01-cause | 8 | 115.213 |

元の 726 固有 ID はすべて保持し、追加した CI 回帰 4 件を内包して 730 件。
7 群の割当と実行 ID の和集合が等しく、重複と欠落が 0 である。
対象の path replay の実成功も JSON に記録した。
C-01 の 3 群は workflow と同じ sparse パターンを適用した checkout で実行した。
Rust と Python の同時実行・5 CPU の競合を含むローカル値であり、Windows 実測へ置き換えない。

C-01 各群の最終値は上表。予備試行を含めた最大観測 157.296 秒を安全側に 2 倍して
Windows の試験部分を約 5 分 15 秒と見積もると、旧 checkout 外れ値 37 分 3 秒、
準備 1 分、CLI build 4 分 32 秒を足して約 48 分。75 分まで約 27 分の余裕がある。
cause-diagnosis 追加 8 件のローカル単独群は予備 151.085 秒、最終
115.213 秒。
Windows での純増は概算約 4〜5 分で、sorted allocation の変化もあるため差分を断定しない。
これは将来 CI の予算推定であり、新先端 Windows は未確定。

PC-3 の固定ソース照合を守るため、既存 runner は元コミットと同一バイトへ戻した。
新規エントリポイントだけが配分を変更する。既存の数値試験・ゲートは変更しない。
実 PowerShell を使う既存試験も実行するため Linux に公式 PowerShell 7.5.4 を配置し、
NC ソース照合・失敗保存の独立 6 件も skip 0 で成功した。
これらと P6-01/P6-03 の追加検査は 730 件に加算しない。

準備中の cargo PATH 不足、fixture の疎配置漏れ、既存 runner 変更を PC-3 が正しく拒否した失敗、
PowerShell 不在で発生した skip は `validation/preliminary/` に保持した。
旧ソース・設定・期待値を変更して通したものではない。再実行を成功件数へ水増ししない。
最終 730 件の各結果と過去試行の区別・正確な件数は `validation.json` を参照する。

新しい validation は、既存の e541fe1 / e8133f0 / fa22e14 と同様に、
先行する不変コミットへ保存した後、P6-01 の追補台帳へ commit、SHA-256、blob SHA-1、
型付き status keys を登録する。過去の登録・証拠を再封印しない。

Windows CIはpush後の新先端では未確定。
C-01は未達、passed=false、原因は判定不能。
