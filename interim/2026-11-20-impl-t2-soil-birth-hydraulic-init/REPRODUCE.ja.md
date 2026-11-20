# T2 の小型施工試験を再現する

PR #13 の当該コミットを、固定診断コミットと証拠アーカイブを含む履歴からチェックアウトする。Windows ではチェックアウト前に `git config --global core.autocrlf false` を設定する。リポジトリ直下で Rust 1.99.1、Python 3.13 と固定依存を使用する。下記の出力先は毎回新規のディレクトリとする。

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$T = 'interim/2026-11-20-impl-t2-soil-birth-hydraulic-init'
$T1 = 'interim/2026-11-20-impl-t1-wet-joint-richards'
$A3 = 'interim/2026-11-26-wet-soil-birth-darcy-a03'
$A4 = 'interim/2026-11-29-richards-wet-joint-a04'
$O = 'runs/t2-new'
python -m pip install -r requirements-verification-win-cp313.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bins --examples
cargo test --locked --release
cargo test --locked --release --examples
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
    python tools/run_ci_python_evidence_shards.py --group $g --output "$O/tests/$g"
}
python "$T/study.py" "$O/controls"
python "$A3/study.py" --steps 4 --output "$O/a03"
python "$A3/water_cache_audit.py" --run "$O/a03"
python "$A3/pointwise_audit.py" --run "$O/a03"
python "$A4/study.py" --steps 4 8 16 --output "$O/a04"
python "$A4/audit.py" --run "$O/a04"
python "$A4/convergence.py" --run "$O/a04" --steps 4 8 16
python "$T1/study.py" a04 "$O/budget-a04" --run "$O/a04"
python "$T1/study.py" refusals "$O/a04-refusals"
python "$A4/joint_controls.py" --output "$O/a04-dry" --kind dry --steps 4
python "$A4/joint_controls.py" --output "$O/a04-hydrostatic" --kind hydrostatic --steps 4
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

Linux は `source /workspace/soil-environment/activate.sh` で固定環境を有効化し、同じ順序で実行する。T2 の `study.py` は施工入力、CLI の標準出力・エラー出力、計算結果 JSON を各ケースの下に保存する。`results.json` の水量は m³、重量は kN、過渡施工の時刻は s。施工時の体積移送、解放途中、最初の正の時間刻みを個別に確認する。詳しい収支行はプロファイル別の `*.json.gz` にある。A-03 の処方ケースには累積水量を作る原票がない。過去の記録は元のコミットで照合する。

Linux 専用のプロセスグループ試験は CI の7群と別に実行する。PowerShell を利用する既存試験には pwsh を PATH に置く。

```bash
python -m pip install -r requirements-reference.lock.txt 'pytest==9.2.0'
python -m pytest tests/linux -q
```
