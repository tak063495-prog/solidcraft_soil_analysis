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
  path replay は switch の index を読むが、switch archive 自体を必要としない。
- 固定コミット・SHA 検証の台帳と他の保存証拠が参照する履歴を守るため、
  `fetch-depth: 0` と blob 遅延取得を維持する。単純な shallow 化は採用しない。
  既存 NC archive の除外パターン、NC の試験・数値設定は変更しない。
- ジョブ上限は **75 分**にする。37 分の checkout の再発も予算に含め、
  完走に必要な余裕を明示的に確保する。cancelled/timeout は集約で引き続き failure になる。
- runner が各試験の実測時間も JSON に記録する。replay の中身・判定・保存値は変更しない。

実際の Git sparse checkout を使う回帰試験で、全群の source 保持と必要 archive の配置、
全試験 ID の一意な割当、CRLF で保存した文書の末尾一致を検査する。

## 検証と時間予算

全件の件数・終了状態・実測時間は `validation.json` と `validation/` の生ログに記録する。
37 分 3 秒の checkout、準備約 1 分、CLI build 最大観測約 4 分 32 秒に、
Windows で実測された mean evidence 約 14 分 16 秒を加えても約 57 分で、
75 分まで約 18 分の余裕がある。Rust integration は同 run 最大約 28 分 12 秒で、
同じ checkout 外れ値と準備を足して約 66 分、余裕約 9 分。
新しい C-01 各群のローカル実測と追加の余裕は最終検証後に追記する。
これらは予算見積りであり、新先端 Windows の完走保証・成功実測ではない。

新しい validation は、既存の e541fe1 / e8133f0 / fa22e14 と同様に、
先行する不変コミットへ保存した後、P6-01 の追補台帳へ commit、SHA-256、blob SHA-1、
型付き status keys を登録する。過去の登録・証拠を再封印しない。

Windows CIはpush後の新先端では未確定。
C-01は未達、passed=false、原因は判定不能。
