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

B-03 の停止理由はリポジトリのルートで個別に確認できる。PowerShell で `$env:PYTHONPATH='tests'` を設定し、次のコマンドを実行する。

```text
python -m unittest -v test_b03_offset_kink_cause_diagnosis.B03OffsetKinkCauseDiagnosis.test_prior_record_and_production_source_bytes_remain_bound
```

停止後の B-03 検査には各診断の起点コミットを含む履歴を用いる。blob の不足時はリポジトリのリモートから取得し、固定コミットに記録された ID と実バイトを照合する。T1 実装前の期待 SHA-256 は B-03 フォルダに保存されている。

```text
git rev-parse 3747630ab615b4be2117ca1aaa13e798cd436374:src/stages.rs
git cat-file -t 8266b160714ec242215048e56d400434f3df5e3e
python -m unittest discover -s tests -p test_b03_offset_kink_cause_diagnosis.py -v
python -m unittest discover -s tests -p test_b03_offset_kink_cloud_followup.py -v
python -m unittest discover -s tests -p test_a03_explicit_shear_initial_stress.py -v
python -m unittest discover -s tests -p test_a03_unequal_thickness_anisotropic_k.py -v
cargo test --locked --release --test wet_joint_richards --test one_way_construction
```

Linux の全件 Python 試験には PC-4 の実プロセス試験を含む。PowerShell のない環境では、既存調査 `interim/2026-12-21-a02-joint-birth-cause-diagnosis/validation/powershell-runtime.json` の Linux x64 版を使用した。配布 tar.gz の SHA-256 は `1fd7983fe56ca9e6233f126925edb24bf6b6b33e356b69996d925c4db94e2fef`。照合後に作業領域へ展開し、pwsh に実行属性を付け、そのディレクトリを PATH に追加する。書込制限のある環境では XDG_CONFIG_HOME、XDG_CACHE_HOME、XDG_DATA_HOME に書込可能な作業用ディレクトリを指定する。pwsh -NoProfile で 7.5.4 を確認してから7群を実行する。初回スキップのログと導入後の実行ログは別々に保存する。

第2回の停止を個別に再現するコマンドを以下に示す。全件コマンドの終了コード101と最初の不一致は evidence/local-verification.zip の rust-all.log と joint-refusal.log に保存した。

```text
cargo test --locked --release --test joint_ext_combos richards_wet_joint_infiltration_drainage_and_two_way_remain_rejected -- --exact --nocapture
cargo test --locked --release -- --list
```

ZIP のバイト列の SHA-256 を evidence/reassessment-at-stop.json の water_ledger_archive、test_logs、acceptance_results_archive と照合してから開く。accepted-water-ledgers.zip の t1、a04、c04 に各基準の JSON 行配列と集計を保存した。accepted-water-ledgers.zip から取り出したケースファイルは gzip 圧縮である。

検証順序の修正後は、旧契約と追加した検査を次のコマンドで別々に再現できる。

```text
git diff 50fb7617 -- src/stages.rs
cargo test --locked --release --test joint_ext_combos richards_wet_joint_infiltration_drainage_and_two_way_remain_rejected -- --exact
cargo test --locked --release --test wet_joint_richards joint_policy_is_checked_before_split_displacement_equilibration -- --exact
cargo test --locked --release --no-fail-fast
cargo test --locked --release --examples --no-fail-fast
```

最終再評価の acceptance.json に記載した source_commit で数値実行を再現する。各群の Python の出力先を分け、全件試験には上記の no-fail-fast を付ける。相対パスと SHA-256 は acceptance.json の archives にある。次のコードはリポジトリルートから保存した3アーカイブの全バイトと ZIP エントリを検査する。

```python
import hashlib, json, zipfile
from pathlib import Path

folder = Path('interim/2026-11-20-impl-t1-wet-joint-richards')
report = json.loads((folder/'acceptance.json').read_text(encoding='utf-8'))
for item in report['archives'].values():
    path = folder/item['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256']
    with zipfile.ZipFile(path) as z:
        assert len(z.infolist()) == item['entries'] and z.testzip() is None
```
