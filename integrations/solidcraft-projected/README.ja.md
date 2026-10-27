# 固定SolidCraftへのMCC画面連携差分（2026-10-27）

固定上流 `edf57cec293dcea301201f1d5173b33e16351466` への62ファイルの全差分です。他の差分と重ねず、クリーンな固定版へ `git apply --check`、続いて `git apply` してください。MIT/Apacheのライセンスと第三者通知を保持しています。上流へのコミットやPRは作成していません。

ConstructionタブでMCC積分法を選び、Save・設計保存・再表示・Undo/Redoを実画面で検証しました。同じ62ファイル版で実CADの8工程を起動し、掘削途中の保存と別ワーカー再開、撤去1134積分点の履歴、連続計算との完全一致、公式GPUの地盤/板/断面/撤去消去を確認しました。生の板値とCAD変更の失効/Undoも検査済みです。

対応する固定ソースで `cargo xtask ci` が合格しました。書式・clippy・全テスト・WASM等を含む結果とソースSHAをmanifestに記録しています。ローカルの実行設定は `SOLIDCRAFT_TEST_JOBS=1`、`RUST_TEST_THREADS=4` です。実画面の実行版も同じ62ファイルに対応します。一般形状・空間収束・全P5/P6・新規PC配布の完了を示すものではありません。

無視されたRust試験は元の設定のままです。Chromeが検出されなかったためWebオフライン試験は省略されており、その範囲も `native-full-ci-scope-and-skips` の報告に明記しています。実デスクトップのGPU確認とは別の記録です。

ビルドと起動時はciプロファイルのパスを明示してください。[属性と起動手順](../../docs/SOLIDCRAFT_ANALYSIS_MODEL.ja.md)、[MCC材料積分](../../docs/MCC_MIDPOINT.ja.md)、[追加検証と残件](../../docs/PROGRESS_20261027.ja.md)を参照してください。
