# ユーザーPCでの再監査・新規求解・残件

対象はPR #13の同じブランチ。PRはマージしない。既存の未保存変更を消さず、作業ツリーを確認して進める。

```sh
git fetch origin codex/plastic-spatial-convergence
git switch codex/plastic-spatial-convergence
git pull --ff-only origin codex/plastic-spatial-convergence
git merge-base --is-ancestor 744fce031c551fe2e1eee9dd69e02f39c2785a9c HEAD
```

最後のコマンドの終了0は指定コミットを含むことの確認であり、数値合格ではない。新規出力先は未存在のパスを使う。保存点のSHA、solver identity、構成則履歴、入力を編集・再封印して再開を通さない。

今回のLinux実行環境はRust 1.99.1、Python 3.13.17、NumPy 2.5.4、SciPy 1.18.1。[environment.json](environment.json) に実測バージョンを記録した。依存関係は既存 `requirements-reference.lock.txt` を使用する。

```sh
rustup toolchain install 1.99.1 --profile minimal --component rustfmt --component clippy
python -m pip install -r requirements-reference.lock.txt
```

Linux/WSLでは `OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 PYTHONUTF8=1 PYTHONIOENCODING=utf-8` をexportする。PowerShellでは次のように設定する。

```powershell
$env:OPENBLAS_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'
$env:CARGO_BUILD_JOBS='2'; $env:RUST_TEST_THREADS='2'
$env:PYTHONUTF8='1'; $env:PYTHONIOENCODING='utf-8'
```

## 分割ファイルと容量

今回の証拠は140,069,972バイト・67分割。再監査だけなら今回の分割ファイルとソース／依存関係でよい。新規求解の追加材料対照は既存oblique証拠282,915,700バイト・135分割を使う。全Python試験には別途、既存mean-stress約5.49 GB、spatial約0.92 GBと、その結合・展開領域が必要である。旧証拠と一時領域には十分な空き容量を確保する。

non-cone疎取得を使って分割ファイルを省いている場合だけ、必要なパスを戻す。

```sh
git sparse-checkout add '/interim/2026-11-10-nc-history-material-cause/*.part*' '/interim/2026-11-10-nc-oblique-and-joint-history/*.part*' '/interim/2026-11-10-nc-mean-stress-release/*.part*' '/interim/2026-11-10-nc-spatial-release/*.part*'
```

今回クラウドではこれらの旧試験・診断用分割を配置して検査した。不足ファイルを理由に試験をスキップしない。通常の全取得ではこのsparseコマンドは不要。新しいcloneでは `--filter=blob:none` を使うと、今回不要な旧巨大証拠までGit objectとして先に取得することを避けられる。

## 1. 生証拠の完全再監査（新しい力学求解なし）

```sh
python interim/2026-11-10-nc-history-material-cause/replay.py --output work/NEW-nc-cause-replay
```

期待出力は `all_parts_and_entries_verified=true`、`parts=67`、`entries=91`、`core_source_files_verified=784`、`sealed_sources_isolated_utf8=true`、`material_history_byte_identical=true`、`whole_saved_runs_exact_except_elapsed_seconds=true`。`numerical_acceptance=false` と `completion_restart_verified=false` は維持される。

全分割・結合ZIP・全エントリーをSHA/サイズ/安全なパスで検査し、封印したソースを `-I -X utf8` の別プロセスで使用する。全773新規受理状態と1,585旧保存場から材料診断を再計測し、全計測JSONのバイト一致を要求する。保存済みの連続・逐次checkpoint・4再開のRunも、計測時間だけを除いて全体一致を要求する。

Windows/別BLASで浮動小数の末尾が変わりバイト比較が失敗した場合は、その失敗を記録する。しきい値を緩めたり、封印した参照を書き換えたりしない。同じ依存関係のLinux/WSLでも再監査できる。これは新しい独立解やWindows上でLinux checkpointを求解再開した証拠ではない。

## 2. 自分のPC・自分の実行版による新規トレースと実再開

```sh
cargo build --locked --release --example nc_history_trace
python interim/2026-11-10-nc-history-material-cause/run_fresh.py --binary target/release/examples/nc_history_trace --output work/NEW-nc-cause-study
```

Windowsでは実行版を `target/release/examples/nc_history_trace.exe` とする。他の引数は同じ。ラッパーは既存の5 kPa物理入力から、撤去初期1/128・1/256の入力を生成する。全受理点のcheckpoint経由トレース、連続実行、除荷0.5／撤去0.5からの4別プロセス再開、旧斜め証拠のSHA照合、全場の材料・接触・力再監査を行う。新しい荷重経路や補助圧力を作らない。

Linuxの今回実測は次のとおり。別PCの最下位桁や停止座標は同一と決めつけず、実出力を確認する。

- `trace128/trace-summary.json`: 受理266、最後の受理λ=0.9512960910797119、切戻し枯渇、最後の真の線形残差1.194127284855861e−10。
- `trace256/trace-summary.json`: 受理507、最後の受理λ=0.9514048099517822、切戻し枯渇、最後の真の線形残差1.0628885795924718e−10。
- `restarts/verification.json`: 2連続対逐次トレースと4別プロセス再開のRun全体が、`elapsed_seconds` だけを除いて一致。
- `measurements.json`: `passed=false`。活動接触0、撤去済み162GP/接触履歴凍結、活動材料の塑性変化、力残差基準1e−8内。全解放・全解放後凍結・独立完了全場・完了再開は未確認。

診断実行版のプロセス終了0は、結果を保存したことを意味する。`complete=false`、cutback、線形ゲート超過を破壊・安定限界・全解放合格へ変更しない。時間予算で止めた場合も別分類として残す。トレース自身の10,000受理点保護で停止した場合は診断エラーであり、破壊や合格ではない。

同一PCで自分の生成した実checkpointを個別に再開する例：

```sh
target/release/examples/nc_history_trace resume work/NEW-nc-cause-study/input128/mesh.json work/NEW-nc-cause-study/input128/ground.json work/NEW-nc-cause-study/trace128/checkpoint-stage-3.json work/NEW-nc-cause-manual-resume
```

Windowsでは先頭の実行ファイル名を `target/release/examples/nc_history_trace.exe` とする。checkpointに記録された実行版・入力と一致しなければ拒否が期待動作である。配布されたLinux checkpointをWindows実行版へ付け替えない。再開Runの比較は `verify_restarts.py` と同じ全体比較を使い、状態や接触項を比較から除外しない。

## 3. 入力のみの独立経路を新しく計算する手順

今回の独立対照は**旧アーカイブの追加診断**であり、新規全解放求解ではない。元の比例撤去を入力だけから新規求解する既存コマンドは次のとおり。

```sh
python tools/low_preload_nc_reference.py --output work/NEW-nc-input-only-128 --steps 128 --method rk4 --budget-seconds 600
python tools/low_preload_nc_reference.py --output work/NEW-nc-input-only-256 --steps 256 --method rk4 --budget-seconds 600
```

`reference.json` の `input_only=true`、`uses_native_fields_or_fractions=false`、`complete`、`failure_reason`、`last_accepted` を確認する。この既存参照は閉鎖固着枝の範囲に限定し、開離・滑りへ出たら拒否する。時間予算・真の線形残差超過等で未完になれば、そのまま未達とする。`passed=false` は既定の安全な分類で、600秒経過を破壊や合格へ変更しない。非比例経路の独立前史400/800・材料内4/8・荷重制御の詳細な既存手順は [bounded-refinement](../2026-11-10-nc-bounded-refinement/REPRODUCE.ja.md) と [oblique](../2026-11-10-nc-oblique-and-joint-history/REPRODUCE.ja.md) を参照する。旧保存場を新しい独立解へ注入しない。

全解放の合格に進めるには、同じ物理入力で元の保持荷重ゼロを確認し、全解放後の実ゼロ増分で全履歴の凍結を検査する。そのうえで、独立塑性3Dの全節点3成分・全材料GP応力/塑性ひずみ/pc/e・全接触履歴、両実装の完了までの刻み感度、実checkpointから完了までの再開一致を揃える必要がある。途中までの変位一致や本体の細分だけでは閉じない。

今回、既存仕様内でこの一式を満たす修正根拠は得られていない。別枝の不存在も証明していないため、材料や接触BIRTH、部分亀裂先端、板／弱結合等の仕様を推測して残件を埋めない。

## 4. Linux/Windowsの全件検証

```sh
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
python tools/run_ci_rust_tests.py --group rust-core --output work/NEW-check/rust-core
python tools/run_ci_rust_tests.py --group rust-integration-0 --output work/NEW-check/rust-integration-0
python tools/run_ci_rust_tests.py --group rust-integration-1 --output work/NEW-check/rust-integration-1
python tools/run_ci_python_tests.py --group python-core-0 --output work/NEW-check/python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/NEW-check/python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/NEW-check/python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/NEW-check/python-spatial-evidence
```

各inventoryの終了0、失敗・エラー・無視・スキップ・Rustフィルター除外0を確認し、Pythonの全発見IDと4群の和集合が欠落・重複なく一致することを確認する。単一群の成功を全件成功としない。Windows固有Rust試験があるためLinuxと件数が異なる。CI成功とNC全撤去の数値受入は別である。
