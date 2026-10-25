# 2026-10-25 途中成果物

P6全体は未完了です。現時点のソース、実斜面の再実行用入力、圧縮した生結果、検証要約を保存しています。

- Rust回帰試験242件、clippy、書式検査が合格しました。
- 固定した224 Tet10・6048積分点の乾燥MCC斜面で、初期平衡→5 kPa載荷→除荷→局所掘削を独立400分割と比較しました。
- 塑性履歴83点、切戻し、3地点の別プロセス再開、非ゼロ体積適合差の保持、材料応力と仕事共役応力の出力を検証しました。
- 掘削底面の掘削増分はRust +4.806507667 mm、独立 +4.806512332 mmです。
- 一般メッシュでの空間収束、他物理・部材・JOINTとの新形式の統合、native画面への新形式接続などが残ります。全要件は `docs/P6_STATUS.ja.md` を参照してください。

## 再実行

リポジトリ直下で、既定の環境構築・依存導入後に実行します。

```powershell
cargo test --release --locked
cargo run --release --locked -- solve-stages interim/2026-10-25/case/mesh.json interim/2026-10-25/case/ground.json interim/2026-10-25/case/case.json result.json
.venv/Scripts/python.exe tools/export_stages.py interim/2026-10-25/case/mesh.json result.json result-vtu
```

`result.json.gz` は当時のRust生結果、`reference-400.json.gz` は独立計算の生結果です。標準gzipで展開できます。当時の実行バイナリSHAは結果の来歴として残し、現在の再ビルド版と同一バイナリだとは扱いません。

`reports/` は公開用の要約と試験ログです。要約には元の記録のSHAを付けています。作業環境の全390証拠ファイル・巨大なビルドディレクトリ・実行環境は同梱していません。既存文書の `runs/` リンクはローカル検証時の来歴です。完成版配布の整備はP6の残件です。
