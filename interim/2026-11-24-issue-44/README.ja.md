# Issue #44　C-01 の再接触点の半滑らか接線

PR #13 の `codex/plastic-spatial-convergence` の先端 `5bd2918814626981053d6880448054c1835a3eaf` を起点とした。対象は `semismooth_contact = true` の Newton 接線と独立 Python の対応する `--active` 経路。[Issue #44](https://github.com/tak063495-prog/solidcraft_soil_analysis/issues/44)

## 原因と修正

開口した確定状態から試行ギャップが厳密にゼロに達すると、直前の接触内分割では閉区間の割合、試行せん断力のノルム、Coulomb 上限が同時にゼロになる。旧実装では `norm > cap` が偽となり、閉側で滑る場合にも固着の微分を返した。Issue の入力での法線列のせん断 x・法線成分は旧式が `(-8, 8)`、実際の閉側が `(-4, 8)`、開側が `(0, 0)`。旧式は二つの片側極限の凸包外にある。

再接触点では、直前の内分割の正のギャップを \(g_0\)、接線増分を \(\Delta u_t\) として、貫入量あたりの試行せん断力の勾配 \(k_t\Delta u_t/g_0\) を求める。大きさが \(\mu k_n\) を超える場合は向きを保って Coulomb 上限の勾配を使う。固着と滑りの勾配が一致する場合も閉側の微分を選ぶ。法線列には閉側の正の剛性を使う。摩擦係数がゼロのときは、閉口側の法線剛性だけを使い、開口側をゼロとする。

既定の Newton 接線は従来の中心差分を使う。接触力・塑性履歴・仕事の更新と出力形式は既存の接触則を使う。受理条件、参照値、既定の工程上限2,000と明示選択の \(B(E)=\lceil2065E/12\rceil\) を保持した。C-01 の全場の刻み収束、漸近域、連続解誤差と他格子での予算検証は未達。[既存の結果と未達](../2026-11-24-c01-event4-retry/README.ja.md)・[旧2,000回の結果](../2026-11-23-c01-method-improvements/README.ja.md)

## 回帰試験

- Rust と独立 Python で、正の摩擦係数の厳密ゼロ再接触について滑りと固着を1・4内分割で検証した。実接触則の開側・閉側それぞれの差分と、選んだ閉側の接線を照合した。
- 滑らかな再接触の固着・滑りを既存の中心差分と照合した。摩擦係数ゼロでは1・4内分割で接線増分がゼロ・微小非ゼロの閉口点、厳密ゼロ再接触、開口点を試し、正の摩擦係数の固着を対照とした。
- Rust の組立試験で、未指定と明示的な `false` の接線を照合し、中心差分と一致することを確認する。能動集合を選択したときも、同じ試行点の状態、節点内力、シリアライズされた接触報告が一致することを確認する。
- 実装修正前は Rust の対象5件中2件が失敗した。独立 Python では滑りの2条件と摩擦係数ゼロの4条件が失敗した。いずれも片側差分または中心差分を照合した後の接線比較で差を検出した。

```bash
cargo test --locked --lib joint_law::semismooth
cargo test --locked --test joints default_recontact_uses_centered_tangent_and_preserves_outputs
python -m pytest -q tests/test_c01_recontact.py
```

## C-01 の再実行

同じ12 Tet10・32接触 GP・7工程で、N1024 の `local-active`、`local-ws`、`ws` を本体と独立 Python で実行する。N512 の `base` を固定上限2,000で再実行し、既定接線の対照と粗変位シードを兼ねる。N1024 は既存の比例上限2,065を使う。各試行の受理履歴と力残差を保存し、既報と工程結果、試行数、終端、全受理の累積塑性散逸を照合する。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
driver=interim/2026-11-23-c01-method-improvements/run_case.py
study=interim/2026-11-24-c01-event4-retry
out=work/issue-44-replay
binary=target/release/examples/c01_methods
cargo build --locked --release --example c01_methods
for solver in native independent; do
  python "$driver" "$out/seeds/$solver-512-base" --solver "$solver" --steps 512 --save-seed --history --binary "$binary"
  for mode in local-active local-ws ws; do
    python "$study/route.py" "$out/routes" --solver "$solver" --steps 1024 --mode "$mode" --seeds "$out/seeds" --binary "$binary" --counts-only
  done
done
```

Windows では例題実行ファイルの拡張子を `.exe` とする。この再実行の比較範囲は工程結果と受理散逸履歴。既報の22場の刻み比較は保存済みの資料を参照。[既報の原票](../2026-11-24-c01-event4-retry/results/routes.json)・[全場の刻み比較](../2026-11-24-c01-event4-retry/results/fields.json)

### 中断時点の再実行結果

完了した5経路は全7工程を完走し、既報と全工程の試行数、結果 JSON、受理散逸履歴の SHA-256 が一致した。

|実装|N|方式|`reload` 試行数|総試行数|
|---|---:|---|---:|---:|
|本体|1024|`local-active`|1,906|9,219|
|独立|1024|`local-active`|1,861|9,170|
|独立|1024|`local-ws`|1,863|9,159|
|独立|1024|`ws`|2,065|9,711|
|独立|512|`base`|1,102|5,056|

独立の粗経路では変位シードの SHA-256 が一致し、過去の粗経路の台帳に記載された結果 JSON の SHA-256 とは差があった。今回の固定上限の粗経路結果 JSON は既報の N512 `base` 本計算原票と一致した。本体の N512 `base` は途中で中断し、そこから作るシードを使う本体の N1024 `local-ws` と `ws` は未実行。

## 検証の中断

`python-core-0` の既存試験2件でエラーを検出し、残りの実行を中断した。部分チェックアウトの準備で、両試験が参照する Git 追跡済みの固定アーカイブ2ファイルを作業ツリーへの配置対象から漏らした。

- `test_cut_volume_history.CutVolumeHistory.test_old_256_path_is_checked_and_read_without_new_solve` が参照する `interim/2026-11-08-wet-joint-cut-volume/evidence.zip`。起点の Git blob は `77a04dda0fe0bc7bc35e8f2566cb67646364dcaa`。
- `test_hydrostatic_cut_volume.HydrostaticCutVolume.test_recorded_256_histories_are_read_with_provenance_not_recomputed` が参照する `interim/2026-11-08-wet-joint-separated-convergence/evidence.zip`。起点の Git blob は `0d60ac27f31cf5cac39501a051ae47ca854adc2e`。

`cargo fmt --all --check` と `cargo clippy --locked --all-targets -- -D warnings` は通過した。Rust の release 試験では追加した接触則4件を含む119件の成功をログで確認した。Rust 全件、既定接線の追加組立試験と doctest は中断時点で未確認。Python は終了済みの実行で、対象と Linux 専用の19件、C-01 専用3群の17件、`python-spatial-evidence` の4件が通過した。Python core 2群と `python-mean-evidence` には完了結果がない。

P6-01 の追補前監査は440ソースで `integrity_ok=true`。今回の固定コミット、台帳追補、共有ブランチへの push、Windows CI、Issue #44 の close は未実施。

### 再開後の確認

指定の cut-volume と separated-convergence の `evidence.zip` を sparse 対象へ追加し、起点の Git blob と一致する内容を `git checkout HEAD` で復元した。前回エラーになった既存試験2件は `python-core-0` の再実行で通過した。

同群の既存試験 `test_t10_reference_tolerance.ReferenceWidths.test_placement_and_both_physical_time_budgets` がエラーになり、実行を中断した。この試験が読む `interim/2026-11-20-impl-t2-soil-birth-hydraulic-init/evidence.zip` も部分チェックアウトの対象から漏れていた。起点の Git blob は `0e6e32c64c307ab781fead794403565399505d4f`。既存の `acceptance.json` に記録された SHA-256 は `e9c3450c32f74cabe79ba1be0db93198efae5c329205763a585c73128febe18a`。

今回の再開分では、後続の Python 群、Rust 全件と既定モードの追加組立試験、最終 clippy、rebase、push、Windows CI、Issue #44 の close には進んでいない。更新後の指示に従い P6-01 台帳への登録は作業対象から外した。

### 追跡ファイル全件の復元

起点 `5bd2918814626981053d6880448054c1835a3eaf` の追跡ファイル13,867件を `git ls-files -t -z` と作業ツリーの実在確認で突き合わせ、未展開の1,035件（1,028個の異なる blob、合計15,491,755,629バイト）を列挙した。`git status` と `git diff --name-status` で既存の変更を特定し、未展開のパスだけを固定コミットから `git checkout --ignore-skip-worktree-bits` で一括復元した。ディスク容量に合わせて15組に分け、各ファイルの `git hash-object --no-filters` を起点ツリーの blob ID と照合した。

全1,035件が一致した。Git の列挙順で「blob ID、タブ、パス、改行」を連結した照合一覧の SHA-256 は `9be7aca468eb00d62220006e2eb7e2383572da5ffdecce31ea3a233d10e15159`。T2 のアーカイブは Git blob と上記の受入記録の SHA-256 の双方が一致した。sparse checkout を解除した後の未展開は0件、ステージ差分は0件。復元した証拠ファイルに変更はない。

### 復元後の検証

`python-core-0` は539件、`python-core-1` は547件を実行し、失敗・エラー・スキップは各0件。前回までにアーカイブ欠落で止まった3件も通過した。Rust の既定モード回帰 `default_recontact_uses_centered_tangent_and_preserves_outputs` は単独実行で1件通過した。

本体の残り3経路も全7工程を完走した。

|実装|N|方式|`reload` 試行数|総試行数|
|---|---:|---|---:|---:|
|本体|512|`base`|1,110|5,077|
|本体|1024|`local-ws`|1,854|9,164|
|本体|1024|`ws`|2,065|9,729|

以上の8経路すべてで、全工程の試行数と終端、結果 JSON と受理散逸履歴の SHA-256 が既報の本計算原票と一致した。N1024 の6経路は入力ファイルの SHA-256 も一致した。N512 は固定上限の粗シード用入力を使っており、本計算原票の比例上限の入力とは SHA-256 が異なる。両実装とも、粗シード票の変位履歴、メッシュ、地盤、施工入力は既報と一致した。粗シード票の結果 JSON には SHA-256 の差があり、本体では既報の粗シード票の `fafc55e2…` に対して今回 `d483c1e8…`、独立では `9e2e6f8b…` に対して `3f6cd970…` だった。いずれも今回の値は対応する N512 本計算原票と一致する。

Python は CI と同じ7群で、重複のない全1,111件を通過した。内訳は core 2群が539件と547件、mean と spatial が各4件、C-01 の switch・path・cause が順に4件・5件・8件。各群とも失敗・エラー・スキップは0件。Linux で mean と spatial の両方が `setUpClass` を完了し、固定アーカイブから再計測した JSON は各原票とバイト単位で一致した。両群の測定記録の `numerical_acceptance` は歴史的な判定どおり `false`。別途 pytest で関連19件（うち CI 群との重複12件、Linux 専用7件）が通過し、Python の実行総数は重複を除くと1,118件。

Rust は `cargo test --locked --release --all-targets` で全540件が通過し、失敗・無視・除外は0件。74個の統合試験ターゲットを含む。既定モード回帰の単独実行1件はこの540件に含まれる。`cargo test --locked --release --doc` は対象0件で正常終了。`cargo fmt --all --check` と `cargo clippy --locked --all-targets -- -D warnings` は通過した。

共有ブランチ上の Windows 検証はコミットごとに実行される。[PR #13 の Checks](https://github.com/tak063495-prog/solidcraft_soil_analysis/pull/13/checks)
