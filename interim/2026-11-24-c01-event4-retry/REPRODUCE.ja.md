# C-01 event4 再試行の再現

12 Tet10 の湿潤7工程をクラウドで実行した。N は空間要素数とは独立した各工程の基準分割数。Python はリポジトリ指定の NumPy/SciPy を含む環境を使用する。作業ディレクトリには新しい空の出力先を指定する。以下はリポジトリ直下から実行する Bash の例。

```bash
set -eu
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
study=interim/2026-11-24-c01-event4-retry
driver=interim/2026-11-23-c01-method-improvements/run_case.py
out=work/c01-event4-retry
seeds="$out/seeds"
runs="$out/runs"
binary=target/release/examples/c01_methods
cargo build --locked --release --example c01_methods
for side in native independent; do
  for n in 128 256 512; do
    python "$driver" "$seeds/$side-$n-base" --solver "$side" --steps "$n" --save-seed --binary "$binary"
  done
  for n in 256 512 1024; do
    for mode in local ws local-ws local-active; do
      python "$study/route.py" "$runs" --solver "$side" --steps "$n" --mode "$mode" --seeds "$seeds" --binary "$binary"
    done
  done
  for n in 512 1024; do
    python "$study/route.py" "$runs" --solver "$side" --steps "$n" --mode base --counts-only --binary "$binary"
  done
  python "$study/route.py" "$runs" --solver "$side" --steps 512 --mode retry --counts-only --binary "$binary"
  python "$driver" "$runs/$side-1024-retry" --solver "$side" --steps 1024 --error-scaled-retry --attempt-budget-policy mesh_proportional --binary "$binary"
done
python "$study/collect.py" "$runs" "$seeds" "$out/results"
```

Windows で個別に実行する場合は実行ファイルに `target/release/examples/c01_methods.exe` を指定する。WS の初期変位には、同一実装、同一物理入力の固定予算 event4 の N/2 解を使用する。今回のクラウド実行では、前回保存した完走解を入力・結果・変位列の SHA-256 と照合して再利用した。再現時も種の生成に要した試行は別経路として扱う。

各経路の `result.json` に工程ごとの試行、受理、棄却理由、最終受理座標と停止理由が入る。棄却した点、未開始工程への補間値は全場出力に含めない。全場は `*-fields/manifest.json` に記載した形状・型の行優先バイナリ。並びは工程順、その中で受理した `i/N` の昇順。集計時にファイルのサイズと SHA-256 を再検査する。刻み提案単独の N1024 試行は工程集計のみを保存した。局所細分を使う3方式では基準点外を含む全受理の32 GP累積塑性散逸も gzip に保存した。符号付き塑性すべりは、3水準共通座標での二つの差ベクトルも little-endian float64 の gzip で同梱した。`fields.json` の `slip_differences` が形状と格納順を示す。全場の生ストリームはこの手順で再生成できる。

```bash
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --release --all-targets --no-fail-fast
cargo test --locked --release --doc
python -m pytest -q tests/test_c01_event4_retry.py
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py --root .
```

Python 全数の対象は `tools/run_ci_python_evidence_shards.py` の7群と `tests/linux`。Windows CI は PR #13 の対象コミットに対して確認する。PC 実行はしていない。PC での再現は上記の12要素・BLAS 1スレッド・1経路ずつとし、前回の計画上限（1経路32 GB・6時間）で打ち切る。全水準の全場ストリームを保存する場合は各経路用のディスク領域も確保する。768 / 6144 は対象外。
