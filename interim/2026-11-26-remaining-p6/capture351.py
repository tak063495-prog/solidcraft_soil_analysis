#!/usr/bin/env python
"""Transparent material-worker wrapper. Capture point 351 for every global evaluation."""
import json
import os
from pathlib import Path
import subprocess
import sys

raw = sys.stdin.buffer.read()
laws, ids, states, strains, counts = json.loads(raw)
i = 351
result = subprocess.run([os.environ["NC_WORKER"]], input=raw, capture_output=True)
row = dict(point=i, law=laws[ids[i]], state=states[i], strain=strains[i], count=counts[i],
           worker_exit=result.returncode, worker_error=result.stderr.decode().strip())
with Path(os.environ["P6_CAPTURE"]).open("a") as stream:
    stream.write(json.dumps(row)+"\n")
sys.stdout.buffer.write(result.stdout)
sys.stderr.buffer.write(result.stderr)
sys.exit(result.returncode)
