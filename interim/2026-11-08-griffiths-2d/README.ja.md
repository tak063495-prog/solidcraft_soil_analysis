# Griffiths斜面661節点比較の記録

2026-11-08。起点は `9187211d4fa224d7eaf936bf73818b06fb9fe190`（PR #13）。新しい斜面求解は実行していない。新h系列の実求解と5％達成、一般3D空間収束、P6は未完了。

元661節点変位は `interim` 内の5つの既存JSONに存在した。固定公開版の元アーカイブを取得し、元の節点座標・拘束・変位との完全一致を確認した。推測・補間は無い。以前の「データが無い」という記述を訂正し、実データを `tests/data/griffiths_2d_reference.json` へ出典付きで保存した。

- `search.json`：起点2639追跡ファイルと101 ZIP／30058メンバーの検索、5つの実データコピー。
- `source-verification.json`：元fixture・ZIP・生結果SHA、全661点の完全一致。Excel再実行・新規求解なし。
- `archived-comparison.json`：既存の実求解変位と元2D実データを3奥行きで再集計。中央面15.581998707／15.885691930／15.593219845／15.672391195％、5％は全て不合格。面内h系列ではなく連続側は判定不能のまま。
- `f1-input-preservation.json`：既定F=1の全3段階の入力4ファイルは指定起点とバイト一致。既存の準備・検査記録も変更なし。
- `matched-input-validation.json`：別系列の参照条件F=1.3500000000000003、1200／4800／19200 Tet10のvalidateのみ。体積220 m³。求解なし。
- `implementation-verification.json`：追加ソース・参照データSHAと、起点の保護対象2101ファイルがバイト一致すること。Rust数値経路と既存失敗記録を変更していない。
- `local-validation.json` と各ログ：fmt、clippy（`--locked --all-targets -- -D warnings`）、release build、全Rust398件、全Python261件（新規22件を含む）が合格。失敗・無視・スキップ0。

比較器は参照ファイルを明示した場合だけ動作する。参照の改変・欠落、F=1との混同、準備manifestや入力の改変、未完成Run、非有限値、物理監査の未達を拒否する。ユーザーPCの新h系列は全点を実節点へ対応させ、既存の独立力・材料・生GP/VTU監査を必須とする。生産CLIの工程メタデータ初期化と、種からの全変位を区別する回帰試験も含む。判定不能と不収束を破壊側上端やFsへ読み替えない。

全Python試験は実装 `7c46773cdb2a0b830c4b4359effcdbe7d47b9632` に対して実行した。Rust本体は起点から不変。本記録の最終追記はログと文書・結果記録だけで、実装を変更していない。

Windows CIは [PR #13の最新checks](https://github.com/tak063495-prog/solidcraft_soil_analysis/pull/13/checks) を参照する。この記録スナップショット作成時点では、最終記録コミットのCI完了前であり、Windows合格を先取りしていない。実行後の判定はPRのchecksと作業報告で対象コミットを確認する。

[日本語のユーザーPC手順・判定範囲](../../docs/GRIFFITHS_2D_COMPARISON.ja.md)。柱のTet4最細6.828％不合格、連続推定1.034％合格、Tet10の5％合格は変更していない。
