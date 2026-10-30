"""Reproduce isolated startup failure and verify actual635 tools imports."""
from pathlib import Path,PurePosixPath
import hashlib,json,os,subprocess,sys,zipfile
p=Path(__file__).resolve().parent
out=p/'fixed34350-import-preflight-v1';out.mkdir();runtime=out/'runtime';runtime.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
quality=p/'Root-source635-promotion-v1/quality-v1';sources=read(quality/'source-snapshot.json')
with zipfile.ZipFile(quality/'source-snapshot.zip') as z:
 assert z.testzip() is None and len(z.namelist())==635 and set(z.namelist())==set(sources)
 for name,h in sources.items():
  q=PurePosixPath(name);assert not q.is_absolute() and '..' not in q.parts
  raw=z.read(name);assert hashlib.sha256(raw).hexdigest()==h
  dest=runtime/name;dest.parent.mkdir(parents=True,exist_ok=True)
  with dest.open('xb') as f:f.write(raw)
probe=runtime/'tools/import-preflight.txt'
probe.write_text("import json\nimport export_stages, verify_stage_package, verify_projected_stage_package, fem_reference, mcc_rate_fem\nfrom pathlib import Path\nmodules=[export_stages,verify_stage_package,verify_projected_stage_package,fem_reference,mcc_rate_fem]\nassert all(Path(m.__file__).resolve().parent==Path(__file__).resolve().parent for m in modules)\nprint(json.dumps({'passed':True,'modules':[m.__name__ for m in modules],'all_modules_from_same_tools_directory':True}))\n",'utf-8')
env=dict(os.environ,PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
for key in ('PYTHONPATH','PYTHONHOME'):env.pop(key,None)
records=[]
for name,flags in [('retained-isolated-failure',['-I','-X','utf8']),('corrected-cli-startup',['-X','utf8'])]:
 command=[sys.executable,*flags,str(probe)]
 result=subprocess.run(command,cwd=runtime,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
  encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW)
 (out/(name+'.log')).write_text(result.stdout,'utf-8')
 records.append(dict(step=name,arguments=command[1:],exit_code=result.returncode))
 if name=='retained-isolated-failure':assert result.returncode!=0 and "No module named 'export_stages'" in result.stdout
 else:assert result.returncode==0 and json.loads(result.stdout)['passed']
assert all(sha(runtime/n)==h for n,h in sources.items())
report=dict(passed=True,actual_Root635_tools_imports_verified=True,records=records,
 same_fixed_Python_environment_same_sources_only_isolation_flag_changed=True,
 original_failure_reproduced_and_retained=True,all635_source_hashes_unchanged=True,
 no_analysis_or_original_data_read_or_written=True,python_probe_not_counted_as_FEM_accuracy_test=True,p6_complete=False)
with (out/'terminal.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
print(json.dumps(dict(passed=True,isolated_failure_reproduced=True,corrected_imports_passed=True)),flush=True)
