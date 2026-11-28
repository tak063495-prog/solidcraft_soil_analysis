PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

所有Popenの終了確認後だけunlockし、確認不能時は所有情報付きlockを保持。通常・起動中signal・待機/status例外と安全な再起動を検証。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`, `1b5fcd86542be615dafcc8dc9b27efca474006b9`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_resume_ownership`：9件。`test_p2_resume_ownership.OwnedChild.test_parent_only_sigint_reaps_child_before_unlock_and_allows_safe_retry` / `test_p2_resume_ownership.OwnedChild.test_parent_only_sigterm_reaps_child_before_unlock_and_allows_safe_retry` / `test_p2_resume_ownership.OwnedChild.test_main_thread_runtime_sigint_reaps_real_child` / `test_p2_resume_ownership.OwnedChild.test_main_thread_runtime_sigterm_reaps_real_child` / `test_p2_resume_ownership.OwnedChild.test_interrupt_during_launch_retains_the_handle_for_cleanup` / `test_p2_resume_ownership.OwnedChild.test_wait_exception_stops_and_waits` / `test_p2_resume_ownership.OwnedChild.test_status_exception_after_launch_stops_and_waits` / `test_p2_resume_ownership.OwnedChild.test_unconfirmed_child_exit_retains_lock_with_owner_record` / `test_p2_resume_ownership.OwnedChild.test_normal_child_failure_unlocks`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[tools/resume_projected_sparse_study.py:128](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/resume_projected_sparse_study.py#L128-L163), [tools/resume_projected_sparse_study.py:197](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/resume_projected_sparse_study.py#L197-L232)。

Linuxの所有子を実行。Windows CIの実子終了回帰はCRT main-thread signalによる経路であり、POSIX外部SIGTERMや実studyの安全停止を認定しない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
