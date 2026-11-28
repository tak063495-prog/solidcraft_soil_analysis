# 正確なFsの証拠検査・実断面SRM準備

起点コミットは `9f6c228d3b40d02b32b2877f0ed3d0feb9279768`。PR #13の既存変更を保持しています。[方法とユーザーPC手順](../../docs/EXACT_FS_REAL_SECTION_SRM.ja.md)。

追加実測は、一様応力・単軸15 kPa・c=10 kPa・φ=ψ=0・6 Tet10（27節点）のみです。epoch/fresh探索とも正確なFs区間は **[1.3320312500000002, 1.3337890625000002]**、幅 **0.0017578124999999556 ≤ 0.002**。完成材料経路の釣合いと平均応力の独立排除を再監査しました。中点をFsにしていません。下端の材料経路進捗は1.0、力残差2.1766882456138767e-11（許容1e-8）、拘束誤差0。上端の平均MC降伏超過0.005125201397193635 kPaは数値余裕1.5008273721818116e-7 kPaを超えています。詳細は `uniaxial-endpoint-evidence.json` です。斜面・実断面のFsではありません。元の[1.3250000000000002, 1.3531250000000001]は変更していません。

`uniaxial-measurement.json` は全実行コマンド・生成solverのSHA・入力/出力のSHA、`uniaxial-evidence.zip` は入力・元Run・全試行結果・判定です。`archived-evidence-summary.json` は公開JSONを読むだけの斜面F=1.45/F=2の要約で、再求解・仮想仕事上端の再計算ではありません。`input-inventory.json` は実断面入力の検索範囲・候補と限界を記録します。

実断面入力の組を特定できず、実断面SRMは未実行・未完了です。準備器は実入力の指定と検査・既存保存ジョブ作成だけを行い、求解や合格を作りません。試験用の小さな一軸ジョブの実行を実断面SRMへ読み替えません。

Linux（Rust 1.99.1、Python 3.13.17、固定9依存）で、fmt、clippy（`--locked --all-targets -- -D warnings`）、release build、全Rust404件・全Python265件が合格しました。失敗・無視は0です。新規試験はRust6件・Python4件。実装コミットは `35d1d52f2e6b2a045abc9e9f1b1a94d452bb90fa`、凍結698ファイルは試験後も一致し、追加測定のsolver SHAは最終release buildと一致しています。詳細は `local-validation.json`、`runtime.json`、各ログ、`local-source-manifest.json`、起点保持は `preservation.json` です。Windows [CI 42418789044](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/42418789044) は成功しました。fmt、clippy、release CLIビルド、Rust406件（Linux404件＋Windows固有2件）、Python265件が合格し、失敗・無視は0です。CIの凍結ソース数は698。実装コミット35d1d52と、CIのマージコミット `10296652219604f2aa93c579e3ec06efe9c0c263` の検証対象ソースはgit比較で差分0です。`windows-verification.json`、`windows-job-status.json`、`windows-test-summary.log` に確認結果を保存しました。artifact13593626257のメタデータを取得しましたが、File Service参照の直接ダウンロードは環境のHTTP 403で取得できず、artifact内部の指紋のバイト比較は未実施です。件数・成功は完了した公開ジョブログと工程状態、ソース対応はGitコミットで確認しています。最終追補はこの証拠ディレクトリのみで、凍結698ファイルは変更していません。

開発中の初回Python4試験では、準備器が既存ジョブmanifestの実ファイル名を `manifest.json` と誤って参照し、計5エラーが発生しました。既存の正しいファイル名 `job.json` へ修正し、新規4試験は通過しています。数値経路には触れていません。初回失敗ログを `retained-preparation-path-failure.log` に保持します。

小ケースのmesh/ground/caseは、既存 `interim/2026-10-30-Root-SRM-and-MCC-recovery/Root-SRM-inputs/closed-form/` の入力をバイト一致でコピーしたものです。fresh SPECの宣言幅0.002は既存 `examples/srm-search.json` と同じです。既存入力と既存区間は変更していません。
