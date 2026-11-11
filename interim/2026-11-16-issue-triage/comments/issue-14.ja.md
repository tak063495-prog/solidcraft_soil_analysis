PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

元のGP非有限値拒否は有効。追記の後続no_flow流量NaNは起点に残存。さらに負/負無限大の相対水収支誤差を実API/CLIで誤受理したため、全境界流量と縮約値の有限scalar検査、相対誤差の有限・非負条件を追加。既存上限は不変。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`, `bfec86f0ff119f44a8046853258e3daf10a992f7`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_soil_source`：2件。`test_p2_soil_source.SoilGP.test_normal_api_and_cli` / `test_p2_soil_source.SoilGP.test_nonfinite_mismatch_and_shape_api_and_cli`
- `test_issue_triage_followups`：3件。`test_issue_triage_followups.SoilBoundaryFlow.test_normal_complete_boundary_inventory_in_all_modes` / `test_issue_triage_followups.SoilBoundaryFlow.test_first_and_last_nonfinite_non_scalar_or_mismatched_flow_are_rejected` / `test_issue_triage_followups.SoilBoundaryFlow.test_reported_balance_error_must_be_finite_scalar_and_nonnegative`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

元のclose理由そのものの回帰ではなく、後続コメントの残存経路を今回解消したものです。修正前の再現と正常対照の差分は記録に残しています。

該当箇所：[tools/verify_soil_source_column.py:71](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/verify_soil_source_column.py#L71-L119)。

実外面12面を使う解析式fixtureと実checker API/CLI。RustがNaNを生成したという指摘ではない。正常reportは修正前後完全一致。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

本issueはclosedを維持し、今回の追記経路の修正根拠を残します。
