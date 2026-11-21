PR #13 の `d8616cefb526e70f422ec8c164afd2353c8f8149` を起点に、同じ `codex/plastic-spatial-convergence` で再照合しました。修正後の判定は **解決済み（本issueの検証不具合の範囲）** です。

旧archive全part/entryのSHA・size・inventoryと可逆JSONLの復元検査を保持し、検証済み旧snapshotを-I -X utf8の別Pythonで再計測。現行checkoutのsource差は診断欄へ分離し、旧実行を新sourceへ結び替えない。

修正コミット：`113d83b8a37b3b68ce46eb00b6664f042334eae0`, `b235824ea55a0e70b8788e586ee8562549f8baf5`, `c0addfb907273a5a71710a2085a487e6e1f90086`。

対応する回帰（全件実行に内包、subTest条件や再実行を重複加算しない）：
- `test_nc_mean_stress_evidence`：4件。`test_nc_mean_stress_evidence.MeanStressEvidence.test_archive_and_exact_raw_json_reconstruction_are_verified` / `test_nc_mean_stress_evidence.MeanStressEvidence.test_exploratory_paths_reach_zero_holding_but_keep_extra_confinement_separate` / `test_nc_mean_stress_evidence.MeanStressEvidence.test_actual_checkpoints_restart_to_completion_with_exact_fields_and_systems` / `test_nc_mean_stress_evidence.MeanStressEvidence.test_forged_acceptance_is_rejected`
- `test_nc_spatial_evidence`：4件。`test_nc_spatial_evidence.SpatialEvidence.test_all_archive_bytes_and_recomputed_measurements_match` / `test_nc_spatial_evidence.SpatialEvidence.test_both_maps_and_orders_have_full_saved_fields_and_frozen_histories` / `test_nc_spatial_evidence.SpatialEvidence.test_separate_process_restarts_match_actual_corresponding_fields` / `test_nc_spatial_evidence.SpatialEvidence.test_forged_acceptance_flag_is_rejected_by_field_audit`
- `test_historical_utf8_isolation`：2件。`test_historical_utf8_isolation.HistoricalUTF8Isolation.test_explicit_utf8_retains_isolation_and_ignores_conflicting_environment` / `test_historical_utf8_isolation.HistoricalUTF8Isolation.test_default_path_read_roundtrips_the_archived_diagnostic_characters`
- `test_spatial_replay_isolation`：3件。`test_spatial_replay_isolation.SpatialReplayIsolation.test_utf8_isolated_child_imports_only_the_selected_snapshot` / `test_spatial_replay_isolation.SpatialReplayIsolation.test_missing_dependencies_are_copied_only_at_the_historical_sha` / `test_spatial_replay_isolation.SpatialReplayIsolation.test_tampered_fallback_is_never_copied`
- `test_nc_replay_contract`：2件。`test_nc_replay_contract.HistoricalReplayContract.test_changed_checkout_uses_verified_snapshot_and_reports_eol_separately` / `test_nc_replay_contract.HistoricalReplayContract.test_corrupt_parts_sources_inventory_and_measurements_are_rejected`

Linuxはfmt・clippy（locked/all-targets、`-D warnings`）・locked release build・全Rust451件・全Python512件（tools含む）が成功し、失敗/無視/skip/除外0です。起点d8616ceのWindows [run 42951479845](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42951479845) は全群成功、Rust453/Python505件。今回push後のWindows成功とは区別します。

同runのmean job 130296554677、spatial job 130296554588で各4試験名の`... ok`と`Ran 4 tests`を確認しました。各1314/411 entriesの検証が実行され、skipによる解決ではありません。Windowsは派生float差を記録し、`measurements_byte_identical=false`、離散証拠一致と受入falseを維持しています。旧part・manifest・source指紋・実行版SHA記録を付け替えていません。

該当箇所：[interim/2026-11-10-nc-mean-stress-release/remeasure.py:21](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/interim/2026-11-10-nc-mean-stress-release/remeasure.py#L21-L25), [interim/2026-11-10-nc-mean-stress-release/remeasure.py:47](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/interim/2026-11-10-nc-mean-stress-release/remeasure.py#L47-L137), [interim/2026-11-10-nc-spatial-release/remeasure.py:23](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/interim/2026-11-10-nc-spatial-release/remeasure.py#L23-L52), [interim/2026-11-10-nc-spatial-release/remeasure.py:73](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/interim/2026-11-10-nc-spatial-release/remeasure.py#L73-L147), [tests/test_nc_mean_stress_evidence.py:15](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/bfec86f0ff119f44a8046853258e3daf10a992f7/tests/test_nc_mean_stress_evidence.py#L15-L55)。

歴史的データの再計測であり旧MCC実行バイナリの再実行ではない。両配布archiveは実行バイナリを省略し、元のSHA/実行記録を保持する。numerical_acceptance/current_checkout_numerical_certificationはfalse。Windowsは派生float差があり計測JSONのバイト一致ではない。

[日本語の全判定・試験・保存照合](https://github.com/tak063495-prog/solidcraft_soil_analysis/tree/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage) / [ToDo](https://github.com/tak063495-prog/solidcraft_soil_analysis/blob/codex/plastic-spatial-convergence/interim/2026-11-16-issue-triage/todo.ja.md)。NC全撤去・湿潤JOINT全履歴収束・曲面Tri6・P6全体は未達のままです。PRはマージしません。

本issueは解決済みとしてcloseします。
