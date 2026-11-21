# 再現手順と試験結果

## 修正前

`repro-before.rs` は本体を変更せずに実行した観測用の Rust 試験。以下は修正前の別作業ツリーで実行する例。

```sh
evidence="$PWD/interim/2026-11-21-issue37-38-wet-joint-transition"
repro_dir="$(mktemp -d)"
git worktree add --detach "$repro_dir/base" 5d0d456393e324685d13964673c9422be1cf5b5a
cp "$evidence/repro-before.rs" "$repro_dir/base/tests/issue37_38_before.rs"
cd "$repro_dir/base"
cargo test --locked --release --test issue37_38_before -- --nocapture
```

実測ログは `repro-before.log`。2 件の `ok` は不具合の再現を確認した結果。#37 は撤去 1 件、16 フレームで完走。#38 は 0.125 s にセル 3・4 の面水頭が 2.00 → 2.25 m に変わり、各セルの全4面と JOINT トレースの外向き流量は前後ともゼロだった。

## 修正後

Rust は `rust-toolchain.toml` の指定を使う。Python は 3.13 と `requirements-reference.lock.txt` の固定依存を使う。試験ランナーは `unittest` と、Linux 専用試験の pytest 9.2.0。

```sh
python -m pip install -r requirements-reference.lock.txt pytest==9.2.0
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --test wet_joint_richards -- --nocapture
cargo test --locked --release --no-fail-fast
cargo build --locked --release
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence \
             python-c01-switch python-c01-path python-c01-cause; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "runs/issue37-38/$group"
done
python -m pytest --version
python -m pytest -q tests/linux/test_pcfail_process_group.py
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

| Python の群 | 割当 | 実行 | 成功 | 失敗 | エラー | スキップ | 未実行 | 終了コード |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `python-core-0` | 485 | 485 | 484 | 0 | 0 | 1 | 0 | 1 |
| `python-core-1` | 498 | 498 | 497 | 0 | 0 | 1 | 0 | 1 |
| `python-mean-evidence` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-spatial-evidence` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-switch` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-path` | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| `python-c01-cause` | 8 | 8 | 8 | 0 | 0 | 0 | 0 | 0 |
| `linux-pytest` | 7 | 7 | 7 | 0 | 0 | 0 | 0 | 0 |

今回の Linux では、mean と spatial の既知のセットアップエラーは再現しなかった。
PowerShell を利用できないため 2 件がスキップされた。記録した終了コードには CI ランナーがスキップを非ゼロとして返す挙動を含む。
unittest の全7群のインベントリーは 1007 件を重複・欠落なしで割り当てた。通常の探索対象外にある Linux 専用試験は pytest 9.2.0 で7件を実行し、Python の総計に加えた。

Rust の最終全件実行は 76 ターゲットの結果を集計し、498 件成功、0 件失敗、0 件無視だった。
過去診断ファイルの展開前に開始した初回の試行は、終了前に中断した。完了したターゲットの小計は成功 215 件、失敗 1 件、無視 0 件で、`rust-initial-incomplete-checkout.log` に残した。失敗した過去診断試験が読む5ファイルを起点コミットから復元し、すべて Git blob SHA を照合した。当該試験の単独再実行と、その後の全件再実行で通過した。

Linux は Rust 1.99.1、Python 3.13.17 を使用した。ロックファイルの Python 依存は全件の版を照合した。
空き容量に合わせて Python の証拠群と core 2群の作業ツリーを順に展開し、各群が参照する証拠を用いて実行した。
P6-01 監査の追補は、既存の追補台帳の全登録を維持し、新しい1件の固定コミット、SHA-256、Git blob SHA-1、型付き状態を検証した。

資料の末尾は UTF-8 の LF で固定した。末尾の文言を検査する Python 試験は、実際に CRLF で書いたファイルを読み直して同じ結果になることも確認する。P6-01 は修正記録の `validation.json` を固定コミット、SHA-256、Git blob SHA-1 で追補台帳に登録する。

FEM の全項目は未達、passed=false。
