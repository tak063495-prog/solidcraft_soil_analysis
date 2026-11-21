# NC の PC 再実行修正と残差・旧証拠の診断

PR #13、ブランチ `codex/plastic-spatial-convergence` の先端 `79d839e3f89331a82d976efc1f23cbe2c848cbee` から作業した。実施日は 2026-11-12（日本時間）。フォルダーの日付は依頼ラベル。PR はマージしない。NC-01 / P6 は未完、数値受入は false のまま。

## 固定入力の比較

`FreshNative` は生成直後の Python 入力と、Rust の型を経て保存された入力を辞書全体で直接比較していた。このため、同じ物理条件でも `gradient` / `gradient_kpa_per_m` / `reference_point_m` のゼロ既定値、空の荷重配列の追加、空 `excavation_sets` と既定 `linear_solver: {method: dense_lu}` の省略を変更として拒否していた。`nc_input_roundtrip` で実際の native スキーマ往復を再現した。

新しい `tools/nc_fixed_inputs.py` はこの既知のスキーマ位置だけに既定値を補い、両側を完全一致で比較する。未知フィールドを削除しない。非ゼロ勾配、全荷重成分、材料、幾何、積分法、Newton 上限、最小刻み、許容値、solver の変更を拒否する。数値比較の許容差は **0**。bool と数値の同一視、NaN / Inf も許さない。比較用コピー以外を変更せず、履歴・受理場・checkpoint・原入力を上書きしない。

`nc_fold_path` の新規診断 checkpoint の source binding に比較器自身の SHA も含めた。既存の荷重全成分・拘束成分・checkpoint 改変・再封印の拒否検査は変更していない。新規 Python 6 件は実 Rust 往復、両刻み、実変更の拒否、高精度残差と記録の保持を検査する。新規 Rust 診断例には非対称系と不正・特異系の 2 件を追加した。

## 実履歴と線形残差

入力生成は既存 midpoint 1/128・1/256 のまま。隔離したソースコピーに、既存 `dense_lu` が拒否した行列・右辺・解を記録する read-only observer を追加した。生産 `src`、方程式、接線、許容値、求解結果は変更していない。全受理点は従来どおり実 `Checkpoint::seal` → JSON → `resume` を通る。両系列の `initial_equilibrated` / `load` / `unload` の全 frame と、撤去前 14 受理点の全行（全 GP・接触履歴を含む）が完全一致した。

|初期刻み|受理点|最終進行度|Linux 最終線形残差|complete / numerical_acceptance|
|---|---:|---:|---:|---|
|1/128|266|0.9514000415802002|1.4875587579753504e-10|false / false|
|1/256|506|0.9514820575714111|1.1446831606412449e-10|false / false|

PC の報告値はそれぞれ線形残差 1.0467e-10 / 1.0273e-10。Linux の残差値まで同一とは主張しない。受理点数と進行度、撤去時のカットバック枯渇を再現した。

拒否された実線形系は 1/128 で 2 件、1/256 で 1 件、すべて 54 自由度。元の行列・右辺は変更せず、native と同じ LU をオフラインで再生した。元の解ベクトルは全 3 件で完全一致した。

|系列・系|native 残差|80 桁 Decimal 残差|2 ノルム条件数|相対後退誤差|反復改良 1 回後の native 残差|
|---|---:|---:|---:|---:|---:|
|128・0|1.2931684994e-10|1.2163225222e-10|5.102684e6|3.727842e-17|4.4382960436e-11|
|128・1（最後）|1.4875587580e-10|1.1816043603e-10|4.486557e6|4.047177e-17|6.3671469403e-11|
|256・0（最後）|1.1446831606e-10|1.0208478989e-10|9.907436e6|3.101773e-17|7.4878730381e-11|

相対後退誤差 `||Ax-b||₂/(||A||F ||x||₂+||b||₂)` は float64 の丸め誤差水準と整合する。一方、高精度評価でも元の解の残差は `1e-10` を超えており、単に残差評価の丸めだけによる誤拒否とは言えない。悪条件化と桁落ちが、右辺で正規化した残差を増幅している。報告する行列積の丸め誤差尺度は約 1.34〜1.42e-8 だが、**受理基準の代用にはしない**。

同じ LU に `b-Ax` の補正を解かせる既存の反復改良は、全 3 系で 1 回目から既存の `1e-10` を満たした。全 8 回の解と残差を保存した。単調減少とは限らない。これは保存系に対する診断であり、本番の Newton や受理済み状態には適用していない。追加の手段は `prepare_nc_linear_observer.py` と `nc_dense_lu_diagnostic` / `diagnose_nc_dense_lu.py` として明示的に実行する。既定求解器の自動補正や許容値の変更は導入していない。これだけで新しい非線形経路・全撤去・fold を受理することはできない。

両方の新しい実履歴は `raw_equal=false`、`fixed_inputs_match=true`、正規化後の差分 0 となった。各履歴から局所 `fold_path` 2 点と `analyze_nc_fold_path` を実行し、すべて rc=0。両者とも `passed=false`、`static_limit_confirmed=false`、`full_removal_reached=false`、点数制限により未確定。128/256 点の長い局所経路を実行したとは数えない。

## 旧 NC mean / spatial の Linux 厳密再計測

A-01 の過去の失敗記録は [そのまま保持](../2026-11-23-joint-combo-extension/validation/legacy-replay-failures/)。既存の Linux バイト完全一致、mean/complement の完全一致、元 measurements、参照値、旧 hash / snapshot、拒否試験は変更していない。Windows の float 差許可は Linux に使用しない。

今回の spatial 全量再計測も `remeasured JSON differs` で不合格。49 個の派生 float、最大 4.440892098500626e-16 の差で、構造・型・離散値の差は 0。再計測 JSON の SHA は `f44ea8b2b5d7420e93b9e14860c51264f052a59d390a6f0cfc4addbcf657d65b` で、A-01 の再計測バイトと完全一致した。元参照の SHA `964215b9ebff7652fa3499607182f097eef07f51a7b1337c117f288e1ca987dc` は維持している。新しい [差分](observations/spatial-replay-differences.json) と再計測 JSON を保存した。

`diagnose_nc_legacy_floats.py` は旧 mean archive の manifest と SHA が一致する controls（`cfe27ea2393ef9c318405f1412f1d730b968a22663c89c4d75f55c1c3b4c9562`）、元 `fem_reference.py`（`79b8602a8a8b9f75d33720db45b383b4ba7476eab515be824950fbe152f937d3`）、関連する元 driver SHA を照合して最小例を再現する。新しい NC 求解や全量再監査の代用ではない。

- 重み 162 成分は完全一致。元 driver と同じ縮約式で再構成した B は 153 成分、最大 4.440892098500626e-16 の差。
- mean/complement の 2 ベクトルは計 10 成分で差があり、最大 2.220446049250313e-16 kN。A-01 の最初の失敗を再現した。
- 元ソースの `Geometry.volume` と standalone `volume` はこの環境で完全一致。キャッシュ導入による演算変更を原因とする証拠はない。C / Fortran / 元状態の 14 列 stride を変えても差は同じ。
- default、Haswell、Sandybridge、Prescott、Nehalem、Zen、および NumPy `X86_V3` 無効化、Prescott との組合せの計 8 条件でも回復しない。

元ソース・入力が同じでも派生 float の厳密一致を再現できず、実行環境の数値演算差と整合する。NC の物理モデルや荷重制御の不具合と断定する根拠は得られていない。元記録は Python / NumPy / SciPy の版を保存しているが、元実行時の BLAS バイナリ・wheel・CPU dispatch の全指紋を保存しておらず、差の発生箇所を一つの環境要因に確定できない。新しい参照値への置換、旧 snapshot の改造、検査の丸め・許容差化をせずに、この環境で元 byte 列へ合わせる修正はできていない。最小診断と全件の厳密失敗を保持し、全件成功へ読み替えない。

## 検証と保存物

**Linux の Python 全件成功は未達。** 旧 NC mean / spatial の setUpClass エラーを維持した上で、最終結果を [validation.json](validation.json) と `validation/` に保存した。

|検証|最終結果|
|---|---|
|cargo fmt --all --check|成功|
|Clippy locked / all-targets / warnings-as-errors|成功|
|Rust locked release all-targets|454 件成功、失敗・ignored・filtered 0（83 集計群）|
|Rust doc test|0 件、成功|
|Python 通常 2 群|276＋284＝560 件成功|
|Python 旧 NC 2 群|各 setUpClass エラー 1、各 4 件の本体試験は未実行|
|Python 全発見 ID|568 件、4 群へ欠落・重複なく割当て、skip / expected failure / unexpected success 0|

Rust の all-targets は従来の 451 件に、従来の example 1 件と新規 example 2 件も含める。初回は大きな証拠の取得・展開と同時実行し、ジョブ復旧の 2 試験が元の 30 秒監視期限へ到達した。`validation/first-tip-rust-all-targets.log` を保持した。I/O 負荷を下げ、`RUST_TEST_THREADS=1` で全件を再実行して成功した。期限・アサーション・数値許容値は変更しない。初回の失敗、単独 smoke、オフライン線形系 3 件、診断経路の点数を最終試験件数へ加算していない。

Python の 8 件未実行を全件成功とは扱わない。参照値や封印 source を変える合格化は採用せず、元の数値 runtime を完全には再現できないという制約と再現手順を残す。現行コードと封印ソースを混同せず、元の生成時 runtime の全情報がないため、元環境由来か元生成時の計測 context 由来かを一意に断定したとは主張しない。

[PC の再実行コマンドと Linux の再現手順](REPRODUCE.ja.md) を参照。完了済み trace を使う場合は節2だけでよく、新規前史は不要。旧 [節3の手順](../2026-11-17-nc-fold-stability/REPRODUCE.ja.md) にも再利用と異常終了検査を追記した。

[観測 manifest](observations/manifest.json)、[集計](observations/summary.json)、`linear-*-systems.jsonl` / `linear-*-diagnosis.json` に元線形系と全反復を保存した。`fold-smoke-evidence.zip` は 2 系列の短経路の全ファイル（材料要求・応答、全場、実線形系、checkpoint と sidecar SHA）と analyze 出力を含む。観測器の元 source SHA・binary SHA も保存した。全 trace の大きな run / accepted / checkpoint は実行 workspace に保持し、ここには入力と各出力の SHA・サイズを記録する。全 trace の可搬アーカイブや新しい受入証拠とは称さない。

案2、正則化、粘性、p′ 下限、成分別解放率、追加拘束、不釣合い力再配分は実装していない。
