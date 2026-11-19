# Issue #36 の Linux 再現と試験

## 修正前と修正後を同じ軽量子で照合

リポジトリルートから Linux の Bash と Python 3.13 で実行する。2 本のランナーと子スクリプトを新しい一時ディレクトリに用意する。

```bash
set -e
base=ad683b30c6cda1f2567a61157d977184731e8fd8
runner=interim/2027-01-08-pc-failure-diagnosis/pcfail_run.py
work=$(mktemp -d)
git show "$base:$runner" > "$work/before.py"
test "$(git hash-object "$work/before.py")" = 1758edc39fe13b10d550c4021cf442ce7d321248
cp "$runner" "$work/after.py"
cat > "$work/child.py" <<'PY'
import subprocess, sys, time
from pathlib import Path
grandchild = subprocess.Popen(
    [sys.executable, "-c", "import time; time.sleep(30)"])
Path(sys.argv[1]).write_text(str(grandchild.pid))
time.sleep(30)
PY
cd "$work"
python - <<'PY'
import ctypes
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

# Adopt the disposable children to reap them even when PID 1 does not.
assert ctypes.CDLL(None, use_errno=True).prctl(36, 1, 0, 0, 0) == 0

def proc_stat(pid):
    try:
        return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
    except FileNotFoundError:
        return None

for version in ("before", "after"):
    for timeout in (10, 1):
        name = f"{version}-{timeout}"
        pidfile = Path(name + ".pid")
        wrapper = subprocess.Popen([
            sys.executable, version + ".py", "--output", name, "--timeout",
            str(timeout), "--", sys.executable, "child.py", str(pidfile),
        ])
        pids = []
        try:
            deadline = time.monotonic() + 5
            while not pidfile.exists() or not pidfile.read_text().strip():
                assert wrapper.poll() is None and time.monotonic() < deadline
                time.sleep(0.01)
            grandchild = int(pidfile.read_text())
            child = int(proc_stat(grandchild)[1])
            pids = [child, grandchild]
            if timeout == 10:
                os.kill(wrapper.pid, signal.SIGINT)
            exit_code = wrapper.wait(timeout=12)
            result = json.loads((Path(name) / "process.json").read_text())
            grandson_stat = proc_stat(grandchild)
            print(name, "exit", exit_code, "status", result["status"],
                  "grandchild", grandson_stat[0] if grandson_stat else "absent")
            assert result["numerical_acceptance"] is False
            assert proc_stat(child) is None
            if timeout == 10:
                assert exit_code == 1 and result["status"] == "failed"
                assert result["error"] == "KeyboardInterrupt()"
                if version == "before":
                    assert grandson_stat[0] == "S"
            else:
                assert exit_code == 2 and result["status"] == "incomplete"
                assert result["timed_out"] and result["exit_code"] == -signal.SIGKILL
                if version == "before":
                    assert grandson_stat is None or grandson_stat[0] == "Z"
            if version == "after":
                assert grandson_stat is None
                try:
                    os.killpg(child, 0)
                except ProcessLookupError:
                    pass
                else:
                    raise AssertionError(f"fixture process group {child} survives")
        finally:
            if wrapper.poll() is None:
                wrapper.kill()
            wrapper.wait()
            for pid in pids:
                try:
                    reaped, _ = os.waitpid(pid, os.WNOHANG)
                except ChildProcessError:
                    continue
                if reaped == 0:
                    os.kill(pid, signal.SIGKILL)
                    os.waitpid(pid, 0)
PY
```

クラウド実測では、修正前の SIGINT はランナー 1108、直接の子 1109、孫 1110 で、終了 1、孫 `S`。タイムアウトはそれぞれ 1112、1113、1114 で、終了 2、孫 `Z`。直接の子は両方で消失し、孫も再現用の親で停止・回収した後に消失した。修正後の SIGINT は 3048、3049、3050、タイムアウトは 3051、3052、3053。終了値はそれぞれ 1 と 2 で、両方とも子・孫の `/proc` エントリが消失し、`killpg(pgid, 0)` は `ProcessLookupError` を返した。

## テストと監査

リポジトリルートで実行する。Python 3.13、`requirements-reference.lock.txt` の依存、pytest 9.2.0、`rust-toolchain.toml` の Rust 1.99.1 を使用した。全件の再現には保存アーカイブが必要。前回の [Linux 実測のパス一覧](../2027-01-09-ci-pcfail-workflow-pin-fix/linux-results.zip) の 1,076 パスを対象 tree の blob ID と照合して配置した。

```bash
python -m pytest --version
python -m pytest -q tests/linux/test_pcfail_process_group.py
cargo build --locked --release
cargo test --locked --release --all-targets
cargo test --locked --release --doc
for group in python-core-0 python-core-1 python-mean-evidence python-spatial-evidence python-c01-switch python-c01-path python-c01-cause
do
  python tools/run_ci_python_evidence_shards.py --group "$group" --output "runs/issue36/$group"
done
python interim/2026-12-06-p6-01-status-reconciliation/audit_status.py
```

各 Python 群の `python-test-inventory.json` で発見 ID、割当、実行件数を照合する。PowerShell ゲート試験は Linux では skip として数え、Windows CI で確認する。skip を含む群で既存 CI ランナーが返す終了 1 は、実行内訳とログを併読する。試験には Linux 専用の pytest 7 件も含める。

## Linux 全件の実測

Rust 全ターゲットは 467 passed、0 failed、0 ignored、0 filtered out。doc-test は 0 件。両コマンドの終了値は 0。

| Python CI 群 | passed | skipped | 終了値 |
|---|---:|---:|---:|
| python-c01-cause | 8 | 0 | 0 |
| python-c01-path | 4 | 0 | 0 |
| python-c01-switch | 4 | 0 | 0 |
| python-core-0 | 476 | 1 | 1 |
| python-core-1 | 481 | 1 | 1 |
| python-mean-evidence | 4 | 0 | 0 |
| python-spatial-evidence | 4 | 0 | 0 |

7 群の合計は 983 固有 ID、981 passed、0 failed・error、2 skipped、0 未実行。skip は `test_nc_pc_gate.RateArcGate.test_powershell_gate_when_available` と `test_nc_pc_gate_schema.PcGateSchema.test_real_rerun_saves_specific_failure_and_preserves_pc3`。Linux 専用 pytest の 7 passed、0 failed、0 skipped を加えた Python 総計は 988 passed、0 failed・error、2 skipped、0 未実行。

台帳の追加行は `register_sources.py` に修正コミット `ea2542403b98c334e957299dcec3dcf4fa965167`、ID `pcfail_issue36_20270110_acceptance`、このフォルダの `acceptance.json` のリポジトリ相対パスを渡して生成した。固定出典を含む 372 出典の監査結果は `integrity_ok=true`、`numerical_acceptance_performed=false`。

A-01〜A-04、B-01〜B-03、C-01〜C-04、NC-01、P6 は未達、passed=false。
