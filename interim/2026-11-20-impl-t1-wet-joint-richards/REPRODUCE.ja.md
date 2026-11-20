# 小型固定格子での再現手順

PR #13 の当該コミットをチェックアウトし、リポジトリのルートから実行する。Rust は rust-toolchain.toml、Python は CI の 3.13 を使用する。Windows ではチェックアウト前に `git config --global core.autocrlf false` を設定する。Python 環境と証拠アーカイブの配置は既存 Windows CI のジョブ定義に従う。以下の出力先には新しい空のパスを選ぶ。各数値経路は12 Tet4を基本とし、既存の乾燥・水没対照も元の小型条件で順番に実行する。

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$T = 'interim/2026-11-20-impl-t1-wet-joint-richards'
$A = 'interim/2026-11-29-richards-wet-joint-a04'
$C = 'interim/2026-11-22-wet-joint-transient-water'
$o = 'runs/t1-new'
python -m pip install -r requirements-verification-win-cp313.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bins --examples
cargo test --locked --release --test wet_joint_richards
python -m unittest discover -s tests -p test_t1_water_balance.py -v
cargo test --locked --release
cargo test --locked --release --examples
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
    python tools/run_ci_python_evidence_shards.py --group $g --output "$o/tests/$g"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
python "$T/study.py" controls "$o/controls"
python "$T/study.py" refusals "$o/a04-refusals"
python "$A/study.py" --output "$o/a04" --steps 4 8 16
python "$A/audit.py" --run "$o/a04"
python "$A/convergence.py" --run "$o/a04" --steps 4 8 16
python "$A/joint_controls.py" --output "$o/dry" --kind dry --steps 4
python "$A/joint_controls.py" --output "$o/hydrostatic" --kind hydrostatic --steps 4
python "$T/study.py" a04 "$o/budget-a04" --run "$o/a04"
python "$C/run_all.py" --output "$o/c04" --probe target/release/examples/wet_joint_transient_probe.exe --solver target/release/soil-fem-3d.exe --mechanical-steps 32 --water-steps 8 16 32 --hydraulic-steps 4 8 16 --timeout 900
python "$C/saturated_control.py" --output "$o/c04" --steps 4 8 16 --timeout 900
python "$C/diagnose.py" --output "$o/c04" --report "$o/c04/cause-diagnostics.json"
python "$T/study.py" c04 "$o/budget-c04" --run "$o/c04"
```

Linux では仮想環境内の python と、同名で拡張子のない実行ファイルを使う。T1 の `study.py` は checkout の examples から入力を作り、CLI の結果から収支を集計する。各ケースの `*-verification.json.gz` と `*-practical.json.gz` は JSON 行配列で、時刻・区分・所有セル／面対／世代と P、N、ΔM、e、D、E、U、ΣD、実際の上下限を記録する。`budget-summary.json` は区分別最大値、最悪の刻み、判定件数。テストのログと Python の発見・割当・実行件数は全件を通して確認する。
