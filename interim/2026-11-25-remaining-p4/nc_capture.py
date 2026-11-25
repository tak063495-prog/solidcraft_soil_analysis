#!/usr/bin/env python
"""Record failing point inputs for the 448-bit unit-shifted P3 worker."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
raw=sys.stdin.buffer.read()
result=subprocess.run([os.environ['NC_WORKER']],input=raw,capture_output=True)
if result.returncode:
    match=re.search(rb'active point (\d+)',result.stderr)
    if match:
        i=int(match[1])
        laws,ids,states,strains,counts=json.loads(raw)
        record=dict(point=i,law=laws[ids[i]],state=states[i],strain=strains[i],count=counts[i],
                    reason=result.stderr.decode().strip())
        with Path(os.environ['NC_FAILURES']).open('a') as f:
            f.write(json.dumps(record)+'\n')
sys.stdout.buffer.write(result.stdout)
sys.stderr.buffer.write(result.stderr)
sys.exit(result.returncode)
