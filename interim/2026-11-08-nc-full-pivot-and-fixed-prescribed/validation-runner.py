import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
sys.path.insert(0, 'tools')
from freeze_sources import freeze, manifest
out=Path('interim/2026-11-08-nc-full-pivot-and-fixed-prescribed')
os.environ.update(PATH='/workspace/.cargo/bin:'+os.environ['PATH'],CARGO_HOME='/workspace/.cargo',RUSTUP_HOME='/workspace/.rustup',CARGO_BUILD_JOBS='2',RUST_TEST_THREADS='2',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
freeze(out)
def run(name,cmd):
 start=time.monotonic()
 with (out/(name+'.log')).open('wb') as f:
  p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 row=dict(name=name,command=cmd,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start)
 text=(out/(name+'.log')).read_text()
 if name=='rust-tests':
  matches=re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',text)
  row['counts']=dict(zip(('passed','failed','ignored','measured','filtered_out'),map(sum,zip(*[list(map(int,m)) for m in matches]))))
 elif name.startswith('python'):
  m=re.search(r'Ran (\d+) tests? in',text);row.update(count=int(m[1]) if m else None,ok=bool(re.search(r'^OK$',text,re.M)),skip_lines=[l for l in text.splitlines() if 'skipped' in l.lower()])
 print(json.dumps(row),flush=True)
 return row
rows=[]
for name,cmd in [('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('build',['cargo','build','--locked','--release'])]:
 rows.append(run(name,cmd))
 if rows[-1]['exit_code']:raise RuntimeError(rows[-1])
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 rust=pool.submit(run,'rust-tests',['cargo','test','--locked','--release','--all-targets'])
 rows.append(run('python-tests',[sys.executable,'-m','unittest','discover','-s','tests','-p','test_*.py','-v']))
 rows.append(run('python-tools-tests',[sys.executable,'-m','unittest','discover','-s','tools','-p','test_*.py','-v']))
 rows.append(rust.result())
source=json.loads((out/'source-snapshot.json').read_text());assert source==manifest()
report=dict(schema='soil-fem-nc-fixed-prescribed-validation/1',commands=rows,all_software_checks_passed=all(x['exit_code']==0 for x in rows),source_file_count=len(source),sources_unchanged_after_tests=True,new_tests=dict(rust=0,python=5),windows_ci='unconfirmed',low_preload_nc_passed=False,p6_complete=False,p3b_p4_rows_complete=False)
(out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
assert report['all_software_checks_passed']
