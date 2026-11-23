# MCC 塑性仕事のクラウド再現

同じチェックアウトのルートでリポジトリ固定の Python 依存を導入する。出力は新しいディレクトリに作る。PC 実行は行っていない。

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
F=interim/2026-11-23-c03-mcc-plastic-work
W=work/new-c03-mcc
cargo build --locked --release --example c03_bulk_cause_probe --example c03_mcc_plastic_work
python "$F/cases.py" "$W" prepare target/release/examples/c03_bulk_cause_probe
for c in m0c1 m1c0 m1c1; do
  python "$F/cases.py" "$W" solve "$c" native || exit 1
  target/release/examples/c03_mcc_plastic_work "$W/basis/mesh.json" "$W/$c/raw" "$W/$c/native-work.jsonl" || exit 1
  python "$F/cases.py" "$W" solve "$c" independent || exit 1
done
python "$F/compare_work.py" "$W" "$W/results"
cargo test --locked --release --example c03_mcc_plastic_work
python -m unittest discover -s tests -p test_c03_mcc_work.py -v
```

比較器は各軸の他入力、工程・進捗の順序、元の力残差、両実装の224時刻を確認する。summary.json には元の実行ファイルの生成した run.json（プロセス・計時情報を含み得る）、工程入力、診断原票のハッシュを記録する。increments.csv の主値・内訳は材料内分割を全て積算した値である。
