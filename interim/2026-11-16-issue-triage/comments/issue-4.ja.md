PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

実際に保存するindent=2の同じバイト列で個別/累積予算を検査。geometryとmanifestも含め、超過時は成功マーカーなし。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_epoch_budget`：2件。`test_p2_epoch_budget.Budgets.test_actual_file_and_total_boundaries_and_no_success_on_overflow` / `test_p2_epoch_budget.Budgets.test_every_payload_uses_same_encoded_bytes_for_size_and_write`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[tools/export_mc_epoch_view.py:17](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/export_mc_epoch_view.py#L17-L31), [tools/export_mc_epoch_view.py:107](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/export_mc_epoch_view.py#L107-L147)。

小型transport fixtureで力学受理部分だけmock。大規模FEM/Native表示の新規受入ではない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
