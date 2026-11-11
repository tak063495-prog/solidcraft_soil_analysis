PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

単調Richardsonを適用する前に連続差分の符号と有限性を確認。正負の振動系列は実CLIでIndeterminate、単調既知次数の対照は維持。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_richardson`：2件。`test_p2_richardson.SignedRichardson.test_oscillating_positive_and_negative_are_indeterminate` / `test_p2_richardson.SignedRichardson.test_monotone_known_orders_and_existing_guards`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[src/plane_strain_2d.rs:569](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/src/plane_strain_2d.rs#L569-L610)。

既存5%と参照値は不変。一般空間収束やGriffithsの未達を完了にしない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
