#!/usr/bin/env python
"""Record the selected point even when the global worker accepts."""
import json
import os
from pathlib import Path
import subprocess
import sys

raw = sys.stdin.buffer.read()
result = subprocess.run([os.environ["NC_WORKER"]], input=raw, capture_output=True)
i = int(os.environ["P7_POINT"])
laws, ids, states, strain, counts = json.loads(raw)
with Path(os.environ["P7_CAPTURE"]).open("a") as stream:
    print(json.dumps(dict(point=i, law=laws[ids[i]], state=states[i], strain=strain[i],
                          count=counts[i], worker_exit=result.returncode,
                          worker_error=result.stderr.decode().strip())), file=stream)
sys.stdout.buffer.write(result.stdout)
sys.stderr.buffer.write(result.stderr)
sys.exit(result.returncode)
