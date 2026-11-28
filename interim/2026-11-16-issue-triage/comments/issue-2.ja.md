PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

構造照査・入力ZIP・exact診断の外部整合条件は最適化で消えない明示例外。normal/-O/PYTHONOPTIMIZE=1で不一致と正常対照を実行。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_optimized_validation`：2件。`test_p2_optimized_validation.OptimizedChecks.test_structural_normal_and_mismatches_never_publish_manifest` / `test_p2_optimized_validation.OptimizedChecks.test_unpack_inventory_bytes_and_hash_in_all_modes`
- `test_p2_acoustic`：1件。`test_p2_acoustic.MaterialInventory.test_exact_changed_input_and_implementation_rejected_in_all_modes`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[tools/beam_reference.py:214](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/beam_reference.py#L214-L239), [tools/solidcraft_structure_fields.py:45](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/solidcraft_structure_fields.py#L45-L78), [tools/mcc_acoustic_exact_report.py:1](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/mcc_acoustic_exact_report.py#L1-L85), [interim/2026-10-27-native-mcc-selector/unpack_inputs.py:1](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/interim/2026-10-27-native-mcc-selector/unpack_inputs.py#L1-L40)。

構造JSON/合成入力・診断CLIの検証。GUIや実CAD/大規模FEMの新しい受入ではない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
