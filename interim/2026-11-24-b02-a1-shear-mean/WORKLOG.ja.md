作業開始時の origin/codex/plastic-spatial-convergence（2026-11-24 JST）　`7a47af617041bf8a409baf77e0602af6cf1fbac0`

- 同先端を基点とする独立した作業ディレクトリで、Tri6 の内部せん断2条件を選択式にした。
- 12×6・18×9 で既定／新方式の2求解を行い、各求解で読み出しを交換した。解析節点場による重心および全積分点のせん断も比較した。
- Linux の `cargo fmt --all -- --check` と `cargo clippy --locked --all-targets -- -D warnings` が成功。T5 曲面殻の Rust 9試験、新方式2試験、既存の4ケース・根元再現2試験、前回の要因調査再現1試験が成功した。
- [A1・A3 の数値と両層の判定](a1.json)・[変更前バイナリとの回帰照合](regression.json)
