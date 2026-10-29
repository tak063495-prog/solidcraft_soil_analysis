# 工程終点指定を含む本体94e70ecのWindows CI

[Windows CI40823809773](https://github.com/tak063495-prog/solidcraft_soil_analysis/actions/runs/40823809773)は全工程が成功しました。artifact12876991725のZIP SHA256は3c61a93af97e1d90040e2fa937aae4b9bc0af2677db6fafe907e9db9cf4c8f05です。

277 Rust試験（ignore0）・171 Python試験、fmt/clippy、実CLIを使用する参照試験が合格しました。配布された全611ソースを94e70ec公開時の凍結一覧と照合し、175ファイルがbyte完全一致、436ファイルはLF/CRLF差だけで、除外0です。Rust1.99.1と固定Python依存も一致しました。この照合でテストを再実行していません。

ci-source-fingerprints.jsonはWindows側、verification.jsonとsource-byte-differences.jsonが全ファイルの照合証拠です。source-snapshot.zip/jsonは今回の文書を含む本体です。94e70ec以後の文書変更と、CIを実行した時点の全611ソースを区別します。全ての非Markdownソースは同じで、数値/CLIコードは変更していません。

[工程終点指定・105456要素MCC](../2026-10-29-stage-endpoints-and-105k/README.ja.md)と、[Native全体CI/降雨の未達記録](../2026-10-29-native-kernel-and-rain-gates/README.ja.md)は、それぞれの条件と版に対応します。今回のFEM CIを、Native65の曲面Boolean失敗の解決や未達精度の受入へ読み替えません。

別CSR候補の2865要素・77355 GP・832板は、同じ元入力の全5工程・全8群一致、表示/力と凍結600ソースによる独立場監査に合格しました。工程終点指定と組み合わせた281 Rust／171 Pythonの品質検査と自身の保存点復旧を別に進めます。まだ本体へ採用していません。

P0〜P6全体は未完です。Native・空間収束・SRM・水理/接触と一般3D/2D照合・新PC/別GPU・完成版配布は継続します。
