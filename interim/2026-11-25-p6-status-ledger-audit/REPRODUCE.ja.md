# 独立監査とLinux試験の再現

リポジトリのルートで実行する。起点は `79d839e3f89331a82d976efc1f23cbe2c848cbee`。
既存ファイルを上書きしない。出力は新しい作業ディレクトリーへ置く。
監査器はPython標準ライブラリとGitだけを使う。既存全件試験には元の固定依存とRust CLI/診断exampleが必要。

```bash
AUDIT=interim/2026-11-25-p6-status-ledger-audit
mkdir -p work/p6-status-check
python "$AUDIT/audit_status.py" --root . > work/p6-status-check/audit-result.json
python "$AUDIT/audit_status.py" --table > work/p6-status-check/ledger.ja.md
python -m unittest discover -s "$AUDIT" -p 'test_*.py' -v
```

監査の終了コードは0=根拠照合errorなし、1=不一致/矛盾/禁止された完了・範囲拡大、2=入力読取エラー。
warning/infoの一覧も確認する。P6合格の判定コマンドではない。
台帳の引用/範囲は人が固定した分類である。新しい証拠が生じたら、元JSONを黙って再生成せず、別のレビュー可能な変更として扱う。

既存Python全件は、既存の収集方法をそのまま使う。新規試験を `tests/` / `tools/` に混ぜない。

```bash
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
export PYTHONDONTWRITEBYTECODE=1
cargo build --locked --release
cargo build --locked --release --examples
python tools/run_ci_python_tests.py --group python-core-0 --output work/p6-status-check/python-core-0
python tools/run_ci_python_tests.py --group python-core-1 --output work/p6-status-check/python-core-1
python tools/run_ci_python_tests.py --group python-mean-evidence --output work/p6-status-check/python-mean-evidence
python tools/run_ci_python_tests.py --group python-spatial-evidence --output work/p6-status-check/python-spatial-evidence
```

既存のCI収集器はtestsとtoolsを全件収集し、4群への一意な割当と実行数を保存する。途中で1群が失敗しても残りの群を実行して記録する。

旧NC mean/spatialの各4件がsetUpClassで止まった場合、通常試験の成功と別に**setupエラー各1、未実行各4**を報告する。
setUpClassエラーは成功でもskipでもない。参照値・整合検査・閾値を変えない。
数値求解を新しく受け入れるためのコマンドではなく、既存ソフトウェア回帰試験を実行する。

全追加差分が指定新規フォルダー内の `A` だけであることを確認する。

```bash
git diff --name-status 79d839e3f89331a82d976efc1f23cbe2c848cbee HEAD
git status --short
```

大規模の数値再実行、実CAD、Richards湿潤JOINT/杭の独立参照作成、実Windows停止、GUI/GPUはこの監査の再現には含めない。
各残件の実行入力・予算・判定は元interimのREPRODUCEとREADMEを参照する。並行draft PRは参照のみ・未マージとして扱う。

今回の生ログ/ZIPは自動承認審査で外部公開が拒否されたため、PRに含めない。公開する集計と失敗理由は `validation.json`、公開ソースの指紋は `source-fingerprints.json`、作業領域に保持する生ログのSHAは `validation-manifest.json` を参照する。上の既存4群のコマンドで新しい作業先へ再現できる。
