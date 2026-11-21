PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

各frameの確定element_material_idsをsampled/exact共通で解決。旧/新Mの候補、定数変更、legacy、非活動を含む不正inventoryを検証。

修正コミット：`5bee3b47f61878e614cae1e7268a57449d2f12d1`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_p2_acoustic`：4件。`test_p2_acoustic.MaterialInventory.test_matset_before_after_and_exact_cli` / `test_p2_acoustic.MaterialInventory.test_material_constants_same_id_and_legacy_fallback` / `test_p2_acoustic.MaterialInventory.test_bad_inventory_rejected_including_inactive_elements` / `test_p2_acoustic.MaterialInventory.test_exact_matset_legacy_and_invalid_inventory_in_all_modes`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

該当箇所：[tools/mcc_acoustic_diagnostic.py:115](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/mcc_acoustic_diagnostic.py#L115-L160), [tools/mcc_acoustic_exact_report.py:30](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tools/mcc_acoustic_exact_report.py#L30-L60)。

診断API/実CLIの合成fixture。実Rust MATSET施工を新たに通したという意味ではない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

既存のclosed判定を維持できます。このファイルは再照合コメントの文案です。
