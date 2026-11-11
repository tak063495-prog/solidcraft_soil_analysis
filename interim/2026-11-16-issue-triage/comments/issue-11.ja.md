PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

全完了recordにも対象ground/caseの未低減工程prefixをverify_originで共通照合。正しく再封印したA/B混在を全phaseと保存検証/復旧で拒否。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `epoch_srm_gate`：2件。`p2_foreign_completed_origin_cannot_borrow_target_identity_at_finish` / `p2_foreign_record_is_rejected_in_all_phases_and_saved_job_recovery`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[src/epoch_srm_gate.rs:402](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/src/epoch_srm_gate.rs#L402-L425), [src/epoch_srm_gate.rs:504](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/src/epoch_srm_gate.rs#L504-L561)。

同一入力の正規保存・復旧一致も既存試験で保持。一般斜面Fsや実断面SRMの受入ではない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
