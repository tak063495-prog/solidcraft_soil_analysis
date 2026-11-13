# 再現・検証手順

起点 `fa22e14701eaa5014b1d8268e768c1ab513db4a6` と今回の保存コミットを比較する。
取消ログの元は run 43446323125 / job 131956600215、比較は run 43427395029 / job 131892091478。
今回も PR #13 のブランチへ push するだけで、マージはしない。

## 環境と全件実行

Rust 1.99.1、Python 3.13.17 と `requirements-reference.lock.txt` の固定依存を用いる。
Linux でも既存 PowerShell 実行試験を省略しないため `pwsh` を PATH に入れる。
今回の Linux 検証では公式 portable PowerShell 7.5.4 を使用した。
PowerShell 不在の skip を全件成功として扱わない。
checkout 前に `git config --global core.autocrlf false` を設定する。
ローカルでは CI で省略していた archive 以外も必要に応じて元の Git blob から配置する。
実際の配置・欠損を確認し、ファイルの省略を試験の skip で回避しない。

```sh
python -m pip install -r requirements-reference.lock.txt
pwsh --version
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release
cargo test --locked --release
cargo test --locked --release -- --list
python tools/run_ci_python_evidence_shards.py --group python-core-0 --output work/validation/python-core-0
python tools/run_ci_python_evidence_shards.py --group python-core-1 --output work/validation/python-core-1
python tools/run_ci_python_evidence_shards.py --group python-mean-evidence --output work/validation/python-mean-evidence
python tools/run_ci_python_evidence_shards.py --group python-spatial-evidence --output work/validation/python-spatial-evidence
python tools/run_ci_python_evidence_shards.py --group python-c01-switch --output work/validation/python-c01-switch
python tools/run_ci_python_evidence_shards.py --group python-c01-path --output work/validation/python-c01-path
python tools/run_ci_python_evidence_shards.py --group python-c01-cause --output work/validation/python-c01-cause
```

7 群の `python-test-inventory.json` の discovered IDs と assigned IDs の和集合を照合し、
重複・欠落が 0、全結果の failure/error/not-run/skip/expected-failure/unexpected-success が 0 であることを確認する。
取消になった path replay は `python-c01-path` にそのまま含まれる。
`tests/test_ci_timeout_configuration.py` は workflow の sparse パターンを一時 Git repository に
適用し、各群の必要 archive と source の実在を検査する。

## P6-01 登録

先行コミットに保存した `validation.json` を変更せず、
`interim/2026-12-10-ci-p601-inventory-fix/register_sources.py` に
40 桁の保存コミットと新規 ID を指定する。
生成された source を既存 `inventory.toml` の末尾へ追記し、
`audit_status.py` の `SUPPLEMENTAL_SHA256` だけを新しい台帳の SHA-256 に更新する。
固定 commit、SHA-256、blob SHA-1、型付き status keys の検査は維持する。

```sh
python interim/2026-12-10-ci-p601-inventory-fix/register_sources.py --commit <保存コミット40桁> --id CI_TIMEOUT_FA22_VALIDATION interim/2026-12-14-ci-timeout-fa22e14/validation.json
python -m unittest discover -s interim/2026-12-06-p6-01-status-reconciliation -p 'test_*.py' -v
python -m unittest discover -s interim/2026-12-07-p6-03-cloud-reconciliation -p 'test_*.py' -v
```

文書の照合は `Path.read_text(encoding='utf-8')` による改行正規化を使う。
CRLF の実バイトを書いた一時ファイルでも同じ末尾が読めることを試験する。
アーカイブと台帳のバイト hash には改行正規化を適用しない。

Windows CIはpush後の新先端では未確定。
C-01は未達、passed=false、原因は判定不能。
