# T10　既存試験の再実行記録

[参照解と許容幅の調査](README.ja.md)の台帳登録コミット `99ed2859f32a273f0bc9876e445fd78797209512` を Linux で検証した。前回ファイル欠測で止まった `test_all_original_rows_and_verified_limited_scopes` は成功。Rust は 512 件、Python は 1049 件成功した。

## 試験前の証拠照合

起点 `1590e85671ed27c2aae82cc18fb4da54d68da71b` から未配置の 474 ファイルを復元した。配置済みの追跡ファイル 13,509 件（31,744,723,108 bytes）を試験対象コミットの Git blob と照合した。リポジトリの Windows CI が全群で除外する旧 NC postfold 分割片 74 件は、この照合の件数に含まない。

対象の T7 `archive-index.json` で全 12 片のサイズ・SHA-256 を照合した。`evidence.zip.part000` の SHA-256 は `3083271d83bc5af9a329ed528d17a7c284c60d6e7bd96b761c31c8529e34fb25`。結合順の全 283,090,557 bytes の SHA-256 は `85b2154a68267a9241d1bc1321b6691ba6f24b1b3233a35c766f1b91e2732b8f`。分割 ZIP として開き、内部 70 件の CRC を検査してから、T8 モジュールの既存 10 件を再実行した。前回の9成功・1エラーと今回の10成功の各ログを収録した。[T7 原票](../2026-11-22-impl-t7-physical-instability/archive-index.json)

## Linux の全件結果

|試験|結果|
|---|---|
|Rust 全ターゲット|104 テストバイナリ、512 成功、失敗・無視・フィルター除外0。doc test 0 件。|
|Python unittest|7群の割当と実行IDの和集合は発見全 1,042 件と一致。重複・未実行0。|
|Linux 専用 pytest|7 成功。|
|P6-01 読取り監査|80 行・412 ソース・171 証拠。`integrity_ok=true`、`numerical_acceptance_performed=false`。|
|書式・静的検査|`cargo fmt --all --check`、`cargo clippy --locked --all-targets -- -D warnings` とも終了0。|

|Python 群|割当|成功|失敗|エラー|スキップ|
|---|---:|---:|---:|---:|---:|
|python-core-0|501|501|0|0|0|
|python-core-1|517|517|0|0|0|
|python-c01-switch|4|4|0|0|0|
|python-c01-path|4|4|0|0|0|
|python-c01-cause|8|8|0|0|0|
|python-mean-evidence|4|4|0|0|0|
|python-spatial-evidence|4|4|0|0|0|

mean と spatial は各4件のメソッドを実行し、両方とも setUpClass エラー0件。T8 の10件、P6-01 の57件、P6-03 の19件は unittest 全件に含む。

## 再現と保存記録

固定版の依存関係と PowerShell 7.5.4 を用い、`cargo build --locked --release --bins --examples` の後、`cargo test --locked --release --all-targets --no-fail-fast` と `cargo test --locked --release --doc` を実行した。Python は `tools/run_ci_python_evidence_shards.py --group <群名> --output <記録先>` を上表の7群すべてで実行し、Linux 専用は `python -m pytest -q tests/linux`。T8 単独の再実行は `python -m unittest discover -s tests -p test_t8_acceptance_criteria.py -v`。環境・実行ID・ログごとのSHA-256はアーカイブ内に収録した。[検証ログと集計](verification.zip)

アーカイブ SHA-256 は `2631ff54f29b47c847e6a554c6eb3a2d0f546c0421492dbe0440ed51fe7a8787`。T10 の固定調査 JSON と本文は当初の登録バイトを保持している。B-01、B-02、B-03、C-04、A-03 は未達。推奨案はユーザー判断待ち。
