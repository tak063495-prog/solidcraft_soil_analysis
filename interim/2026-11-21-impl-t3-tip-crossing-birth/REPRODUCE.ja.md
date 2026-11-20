# T3 の再現手順

PR #13 の履歴と A-02 の分割アーカイブを含むチェックアウトを使う。Rust は rust-toolchain.toml の 1.99.1、Python は 3.13。Windows ではチェックアウト前に core.autocrlf を false に設定する。作業ディレクトリはリポジトリ直下、結果の出力先は毎回新規にする。

```powershell
$env:CARGO_BUILD_JOBS = '2'
$env:RUST_TEST_THREADS = '2'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
$T = 'interim/2026-11-21-impl-t3-tip-crossing-birth'
$O = 'runs/t3-reproduction'
python -m pip install -r requirements-verification-win-cp313.lock.txt
cargo fmt --all --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release --bins --examples
cargo test --locked --release --test t3_tip_birth
python "$T/study.py" "$O/cases"
cargo test --locked --release
cargo test --locked --release --examples
foreach ($g in 'python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause') {
    python tools/run_ci_python_evidence_shards.py --group $g --output "$O/tests/$g"
}
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
python -m unittest discover -s interim/2026-12-06-p6-01-status-reconciliation -p test_status_audit.py
```

Linux は requirements-reference.lock.txt と pytest==9.2.0 を pip で導入し、同じ7群と `python -m pytest -q tests/linux` を実行する。PowerShell を使う既存試験には pwsh 7.5.4 を PATH に追加する。ホームに書き込めない環境では XDG_CONFIG_HOME、XDG_CACHE_HOME、XDG_DATA_HOME を作業領域に設定する。

study.py は 11-25 A-02 の evidence-manifest.json に従い、各 part、結合 ZIP、取り出す入力メンバーの SHA-256 を照合する。元の mesh、ground、case、投影要求をそのまま本体 CLI と公開 API 診断に渡す。new.log の終了コード1は、拒否理由、結果ファイルの有無と合わせて読む。入力不足や診断器の異常終了も終了コードだけで合否にしない。

S1・S2 の独立側は入力から体積 GL125、両面の Tri6 積分、拘束零空間を組み立てる。初期平衡の終点を直接解き、以後は入力の各工程を4等分して、配置、撤去、摩擦、節点力を計算する。本体の解・受理列・接線は独立求解への入力にしない。比較時に本体の工程末フレームと途中チェックポイントを読み、工程の中点からは別プロセスで最後まで再開する。全 Run の比較から除くキーは elapsed_seconds。出力の results.json に各時点の生差、凍結・新設、拒否診断、実行版 SHA-256 がある。

各施工経路は順番に実行し、既存の81節点・24 Tet10・648材料積分点・32接触積分点に限定する。大規模な細分化は再現に必要ない。
