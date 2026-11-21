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

Rust は `rust-toolchain.toml` の指定を使う。Python は 3.13 と `requirements-reference.lock.txt` の固定依存を使う。試験ランナーは `unittest`。

```sh
python -m pip install -r requirements-reference.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --test wet_joint_richards -- --nocapture
cargo test --locked --release --no-fail-fast
cargo build --locked --release
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence \
             python-c01-switch python-c01-path python-c01-cause; do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "runs/issue37-38/$group"
done
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

<!-- linux-details -->

資料の末尾は UTF-8 の LF で固定した。末尾の文言を検査する Python 試験は、実際に CRLF で書いたファイルを読み直して同じ結果になることも確認する。P6-01 は修正記録の `validation.json` を固定コミット、SHA-256、Git blob SHA-1 で追補台帳に登録する。

FEM の全項目は未達、passed=false。
