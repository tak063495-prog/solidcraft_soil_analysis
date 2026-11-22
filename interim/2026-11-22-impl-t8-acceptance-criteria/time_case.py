"""Record Linux/WSL elapsed time, child peak RSS (KiB) and exit status for one route."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import time


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("report", type=Path)
    p.add_argument("command", nargs=argparse.REMAINDER)
    args = p.parse_args()
    with args.report.open("x", encoding="utf-8", newline="\n") as f:
        start = time.monotonic()
        try:
            result = subprocess.run(args.command, timeout=86400)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
        json.dump(dict(seconds=time.monotonic()-start,
                       child_maxrss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                       exit_code=code), f, indent=2)
        f.write("\n")
    raise SystemExit(code)
