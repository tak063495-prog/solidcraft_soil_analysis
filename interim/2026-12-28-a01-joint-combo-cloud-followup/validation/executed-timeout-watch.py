import os,signal,time,json
from pathlib import Path
pid=23857
cmd=Path('/proc')/str(pid)/'cmdline'
if cmd.exists():
 text=cmd.read_bytes().decode().replace('\0',' ')
 assert 'a01follow_run.py' in text and 'stacked__independent__Pindependent__Rnative' in text
 stat=(Path('/proc')/str(pid)/'stat').read_text().split();start=float(stat[21])/os.sysconf(os.sysconf_names['SC_CLK_TCK'])
 deadline=start+900
 while cmd.exists() and time.monotonic()<deadline:time.sleep(min(30,max(0.1,deadline-time.monotonic())))
 if cmd.exists():
  os.kill(pid,signal.SIGTERM)
  out=Path('work/a01follow-stacked-extra/stacked__independent__Pindependent__Rnative__g8')
  result=dict(complete=False,diagnostic_only=True,pid=pid,predictor='independent',repartition='native',accepted_states=sum(1 for _ in (out/'fields.jsonl').open()),elapsed_seconds=time.monotonic()-start,
   failure_reason='900-second diagnostic wall-clock budget exhausted; partials retained',numerical_outcome='incomplete',scope_refusal=False,physical_failure=False)
  (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
