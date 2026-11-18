# 再現手順

リポジトリ直下から Git 2.52、プロジェクトの Rust toolchain と固定 Python 依存を使う。Windows は checkout 前に `git config --global core.autocrlf false` を設定する。Git のパス引数はリポジトリ相対の POSIX 表記を使う。SHA256 は実バイト、文書の末尾は UTF-8 の universal newline で検査する。

```sh
git diff 18d54b1^ 18d54b1 -- .github/workflows/fem-ci.yml
git diff 18d54b1 20a6940 -- .github/workflows/fem-ci.yml
cargo build --locked --release --all-targets
cargo test --locked --release --all-targets
cargo test --locked --release --doc
python -m unittest discover -s tests -p 'test_*.py' -v
python -m unittest discover -s tools -p 'test_*.py' -v
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

Python は tests と tools の両方が対象。PowerShell ゲートは Linux で skip、Windows CI で確認する。再現時は対象コミットの保存アーカイブを配置し、取得した各ファイルの `git hash-object` と `git rev-parse <commit>:<path>` の出力を照合する。P6 登録は先に acceptance.json をコミットし、register_sources.py にその 40 桁 SHA と相対パスを渡す。出力を追補台帳にバイトで追記し、台帳全体の SHA256 を監査器に固定する。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
