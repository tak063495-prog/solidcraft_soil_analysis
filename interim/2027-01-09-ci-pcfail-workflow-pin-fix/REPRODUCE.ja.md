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

## この実行で使用した取得経路と試験群

次の例は同じ対象を tree、Git 転送、GitHub の実ファイルで照合する。ファイルの配置先には新しい作業ディレクトリを使う。

```sh
rev=20a69407bee17bdeebebdb7e10025f9e91a9f9f3
path=interim/2026-11-10-nc-spatial-release/evidence.tar.gz.part054
oid=$(git rev-parse "$rev:$path")
git rev-parse --is-shallow-repository
git check-attr filter -- "$path"
git fetch origin "$rev"
git fetch --refetch --no-filter --no-tags origin "$oid"
curl -fLsS "https://raw.githubusercontent.com/tak063495-prog/solidcraft_soil_analysis/$rev/$path" -o sample.blob
test "$(git hash-object sample.blob)" = "$oid"
```

Linux の通常集合は CI と同じモジュール単位の 7 群で実行した。各群の `python-test-inventory.json` に発見・割当・実行 ID、失敗・error・skip 件数が出る。既存の CI ランナーは skip があれば終了 1 を返すため、件数と個別ログを併読する。

```sh
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause
do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "runs/pin-fix/$group"
done
```

集計で照合する実測は Rust 467 成功、Python 979 成功・2 skip。両言語とも失敗・アーカイブ不足の未実行は 0。validate 用の固定出典は同じ登録コマンドで validation.json の保存コミットから生成する。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
