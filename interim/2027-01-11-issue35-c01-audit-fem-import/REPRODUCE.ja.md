# Issue #35 の再現・検証

Python 3.13 と `requirements-reference.lock.txt` の依存、pytest 9.2.0、Rust 1.99.1 を使用した。コマンドはリポジトリルートから実行する。

```text
python -m pip install -r requirements-reference.lock.txt "pytest==9.2.0"
```

## 修正前の実測

Bash で修正前のソースを一時 worktree に展開し、修正コミットの回帰試験から正常入力の 2 件を実行する。一時 worktree には Python ソースだけを配置する。

```bash
base=8cf97d2bb4c3e5b25667445f6c0033379d565ec1
fix=04bc2e6d128289febd82fdd2e5a87f5103a83d11
audit=interim/2026-12-11-c01-global-step-control/audit.py
test_path=tests/test_c01_global_field_audit.py
before="$(mktemp -d)/issue 35 before"
git worktree add --detach --no-checkout "$before" "$base"
git -C "$before" sparse-checkout set --no-cone \
  '/interim/**/*.py' '/tools/*.py' '/tests/test_c01_global_step_control.py'
git -C "$before" reset --hard HEAD
git -C "$before" rev-parse "HEAD:$audit"
git show "$fix:$test_path" > "$before/$test_path"
(cd "$before" && python -m pytest -q "$test_path" -k minimal --tb=short)
```

監査ソースの blob は `c270508abe7a4e9f8286cd3c2320a4f632315333`。両ケースとも `study.py:50` のラッパー呼出しから `audit.py:38` の `fem.elastic_tensor(...)` に達した。2 failed、0 skipped、3 deselected、終了値 1。

```text
E   NameError: name 'fem' is not defined
```

## 修正後の回帰試験

次のコマンドは Linux と Windows PowerShell で実行できる。合成データは Tet4 体積要素 1 個、接触積分点 1 個、各工程 1 状態、合計 7 状態。native と independent の履歴を同じ `Writer` で作る。釣合い演算子を合成し、ソルバーを起動せずに実際の `study.audit` を実行する。

```text
python -m pytest -q "tests/test_c01_global_field_audit.py"
```

正常な両経路の最大再組立力残差と最大弾性応力整合誤差は 0。従来の力残差許容値は native が `1e-8`、independent が `1e-10`。応力 1 成分を 1 kPa だけずらすと、既存の応力照査で `AssertionError` が発生する。旧ソースの SHA-256 は 2 本の保存指紋と固定コミットのバイト列で照合し、現行ソースの Git blob ID と SHA-256 は今回の記録と照合する。文書末尾は LF と CRLF の両方で検証する。 [修正前後の指紋](acceptance.json)

## Linux 全件と P6-01

全件試験では過去の Linux 実測で使った 1,076 パスを起点 tree の blob ID と照合して配置する。PowerShell ゲート試験の skip は個別に計上する。 [保存アーカイブのパス一覧](../2027-01-09-ci-pcfail-workflow-pin-fix/linux-results.zip)

```bash
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
export CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2
export PYTHONUTF8=1 PYTHONIOENCODING=utf-8
cargo test --locked --release --all-targets
cargo test --locked --release --doc
cargo build --locked --release
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause
do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "runs/issue35/$group"
done
python -m pytest -q "tests/linux/test_pcfail_process_group.py"
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

P6-01 の追補用エントリは既存のスクリプトから生成した。旧台帳に追記し、`audit_status.py` の追補台帳 SHA-256 を `85e69d9683c6d1d70317795c8d1fe622d269bd2df43ef0b3a142245ddd6fb6a9` に更新した。

```bash
python interim/2026-12-10-ci-p601-inventory-fix/register_sources.py \
  --commit 04bc2e6d128289febd82fdd2e5a87f5103a83d11 \
  --id c01_global_step_issue35_20270111_acceptance \
  interim/2027-01-11-issue35-c01-audit-fem-import/acceptance.json
```

### 全件の内訳

Rust 全ターゲットは 467 passed、0 failed、0 ignored、0 filtered out。doc-test は 0 件。両コマンドの終了値は 0。

| Python CI 群 | passed | failed・error | skipped | 終了値 |
|---|---:|---:|---:|---:|
| python-core-0 | 468 | 0 | 1 | 1 |
| python-core-1 | 494 | 0 | 1 | 1 |
| python-mean-evidence | 4 | 0 | 0 | 0 |
| python-spatial-evidence | 4 | 0 | 0 | 0 |
| python-c01-switch | 4 | 0 | 0 | 0 |
| python-c01-path | 4 | 0 | 0 | 0 |
| python-c01-cause | 8 | 0 | 0 | 0 |

7 群の合計は 988 固有 ID、986 passed、0 failed・error、2 skipped、0 未実行。既存の CI ランナーは skip を含む群で終了値 1 を返す。skip は以下の 2 件。

- `test_nc_pc_gate.RateArcGate.test_powershell_gate_when_available`
- `test_nc_pc_gate_schema.PcGateSchema.test_real_rerun_saves_specific_failure_and_preserves_pc3`

Linux 専用 pytest の 7 passed、0 failed、0 skipped を加えた Python 合計は 993 passed、0 failed・error、2 skipped。今回の回帰試験だけを pytest で実行すると 5 passed、0 failed、0 skipped、6 サブテスト成功。Python CI 7 群との重複分は合計に加えない。

台帳監査は 373 出典で `integrity_ok=true`、`numerical_acceptance_performed=false`。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
