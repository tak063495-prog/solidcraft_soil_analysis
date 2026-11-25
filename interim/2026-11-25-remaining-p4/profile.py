"""Enforce the per-route wall-time cap and record Linux child peak RSS."""
import json
import resource
import subprocess
import sys
import time

start = time.monotonic()
result = subprocess.run(['timeout','24h',*sys.argv[1:]])
print(json.dumps(dict(exit_code=result.returncode,seconds=time.monotonic()-start,
                      peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)))
sys.exit(result.returncode)
