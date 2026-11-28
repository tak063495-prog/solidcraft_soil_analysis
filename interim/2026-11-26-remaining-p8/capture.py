#!/usr/bin/env python
"""Record both requested material inputs from the same assembled global evaluation."""
import json
import os
import subprocess
import sys

raw = sys.stdin.buffer.read()
proc = subprocess.run([os.environ["NC_WORKER"]], input=raw, capture_output=True)
laws, ids, states, strain, counts = json.loads(raw)
with open(os.environ["P8_CAPTURE"], "a") as stream:
    records = [dict(point=i, law=laws[ids[i]], state=states[i], strain=strain[i],
                    count=counts[i]) for i in (351,272)]
    print(json.dumps(dict(exit=proc.returncode, error=proc.stderr.decode(), points=records)),file=stream)
sys.stderr.buffer.write(proc.stderr)
sys.stdout.buffer.write(proc.stdout)
sys.exit(proc.returncode)
