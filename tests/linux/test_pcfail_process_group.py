"""Exercise the diagnosis runner with disposable Python processes on Linux."""
import ctypes
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest


RUNNER = (
    Path(__file__).resolve().parents[2]
    / "interim/2027-01-08-pc-failure-diagnosis/pcfail_run.py"
)
CHILD = """
import os
from pathlib import Path
import subprocess
import sys
import time

pidfile = Path(sys.argv[1])
mode = sys.argv[2]
grandchild = subprocess.Popen([
    sys.executable, "-c", "import time; time.sleep(30)" if mode in
    {"waiting", "orphan"} else "pass",
])
temporary = pidfile.with_suffix(".tmp")
temporary.write_text(f"{os.getpid()} {grandchild.pid}", encoding="ascii")
temporary.replace(pidfile)
if mode == "waiting":
    time.sleep(30)
elif mode != "orphan":
    sys.exit(int(mode) if grandchild.wait() == 0 else 99)
"""
WAIT_FAILURE = """
import os
from pathlib import Path
import runpy
import sys
import time

sys.argv.pop(0)
pidfile = Path(sys.argv.pop(0))
original_waitid = os.waitid
def fail_after_spawn(*args):
    os.waitid = original_waitid
    deadline = time.monotonic() + 5
    while not pidfile.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    if not pidfile.exists():
        raise AssertionError("child never became ready")
    raise RuntimeError("injected wait failure")

os.waitid = fail_after_spawn
runpy.run_path(sys.argv[0], run_name="__main__")
"""


@unittest.skipUnless(sys.platform == "linux", "requires Linux process groups and /proc")
class PcfailProcessGroup(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="pcfail process group ")
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.pidfile = self.work / "pids"
        self.child = self.work / "child.py"
        self.child.write_text(CHILD, encoding="utf-8")
        self.output = self.work / "result"

        # Adopt any fixture orphan left by a failing runner in this test only.
        libc = ctypes.CDLL(None, use_errno=True)
        previous = ctypes.c_int()
        self.assertEqual(libc.prctl(37, ctypes.byref(previous), 0, 0, 0), 0)
        self.assertEqual(libc.prctl(36, 1, 0, 0, 0), 0)
        self.addCleanup(lambda: self.assertEqual(libc.prctl(36, previous.value, 0, 0, 0), 0))

    def start(self, mode, *, timeout=10, fail_wait=False):
        command = [sys.executable]
        if fail_wait:
            command += ["-c", WAIT_FAILURE, str(self.pidfile)]
        command += [
            str(RUNNER), "--output", str(self.output), "--timeout", str(timeout),
            "--", sys.executable, str(self.child), str(self.pidfile), mode,
        ]
        with (self.work / "wrapper.log").open("w", encoding="utf-8") as log:
            wrapper = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        self.addCleanup(self.clean_fixture, wrapper)
        deadline = time.monotonic() + 5
        while not self.pidfile.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        self.assertTrue(self.pidfile.exists(), (self.work / "wrapper.log").read_text())
        return wrapper, self.fixture_pids()

    def fixture_pids(self):
        return [int(pid) for pid in self.pidfile.read_text(encoding="ascii").split()]

    def clean_fixture(self, wrapper):
        if wrapper.poll() is None:
            wrapper.kill()
        wrapper.wait(timeout=5)
        if not self.pidfile.exists():
            return
        for pid in self.fixture_pids():
            try:
                reaped, _ = os.waitpid(pid, os.WNOHANG)
            except ChildProcessError:
                continue
            # A waitable child cannot have its PID reused while we own its reap.
            if reaped == 0:
                os.kill(pid, signal.SIGKILL)
                os.waitpid(pid, 0)

    def check_result(self, wrapper, pids, *, cli_exit, status, **fields):
        wrapper_exit = wrapper.wait(timeout=8)
        self.assertEqual(wrapper_exit, cli_exit, (self.work / "wrapper.log").read_text())
        result = json.loads((self.output / "process.json").read_text(encoding="utf-8"))
        self.assertEqual(result["status"], status)
        self.assertIs(result["numerical_acceptance"], False)
        for key, value in fields.items():
            self.assertEqual(result[key], value)
        for pid in pids:
            self.assertFalse(Path(f"/proc/{pid}").exists(), f"fixture process {pid} remains")
        with self.assertRaises(ProcessLookupError):
            os.killpg(pids[0], 0)

    def test_sigint_to_runner_pid_only(self):
        wrapper, pids = self.start("waiting")
        wrapper.send_signal(signal.SIGINT)
        self.check_result(wrapper, pids, cli_exit=1, status="failed", error="KeyboardInterrupt()")

    def test_exception_while_waiting(self):
        wrapper, pids = self.start("waiting", fail_wait=True)
        self.check_result(
            wrapper, pids, cli_exit=1, status="failed", error="RuntimeError('injected wait failure')",
        )

    def test_timeout(self):
        wrapper, pids = self.start("waiting", timeout=1)
        self.check_result(
            wrapper, pids, cli_exit=2, status="incomplete", exit_code=-signal.SIGKILL, timed_out=True,
        )

    def test_completed_child(self):
        wrapper, pids = self.start("0")
        self.check_result(wrapper, pids, cli_exit=0, status="done", exit_code=0, timed_out=False)

    def test_nonzero_child(self):
        wrapper, pids = self.start("7")
        self.check_result(wrapper, pids, cli_exit=1, status="failed", exit_code=7, timed_out=False)

    def test_incomplete_child(self):
        wrapper, pids = self.start("2")
        self.check_result(wrapper, pids, cli_exit=2, status="incomplete", exit_code=2, timed_out=False)

    def test_child_exits_first(self):
        wrapper, pids = self.start("orphan")
        self.check_result(wrapper, pids, cli_exit=0, status="done", exit_code=0, timed_out=False)
