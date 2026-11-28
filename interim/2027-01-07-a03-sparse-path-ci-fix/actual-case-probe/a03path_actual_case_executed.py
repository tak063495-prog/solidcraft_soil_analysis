"""Run actual old/fixed sparse tests with Windows relative paths on Linux."""
from pathlib import Path, PureWindowsPath
from unittest.mock import patch
import hashlib,json,subprocess,types,unittest
root=Path.cwd();name='tests/test_a03_unequal_thickness_anisotropic_k.py';target=root/name
out=root/'interim/2027-01-07-a03-sparse-path-ci-fix/actual-case-probe';out.mkdir()
source_type=type(root);original=source_type.relative_to;rows=[]
for version,source in [('old',subprocess.check_output(['git','show','6ee03451c49b86bc9f238ba0bed2f21b15f4299c:'+name])),('fixed',target.read_bytes())]:
    namespace=types.ModuleType('a03path_'+version);namespace.__file__=str(target)
    exec(compile(source,str(target),'exec'),namespace.__dict__)
    seen=[]
    def relative(path,*args,**kwargs):
        value=original(path,*args,**kwargs)
        if value.parts and value.parts[0]=='interim' and '.part' in value.name:
            value=PureWindowsPath(value.as_posix());seen.append(str(value))
        return value
    result=unittest.TestResult()
    with patch.object(source_type,'relative_to',relative):
        namespace.A03UnequalAnisotropic('test_ci_sparse_selects_new_archive').run(result)
    assert result.testsRun==1 and not result.errors and not result.skipped
    assert len(result.failures)==(1 if version=='old' else 0)
    assert len(seen)==(10 if version=='old' else 20)
    text=''.join(f[1] for f in result.failures) or 'OK: exact archive byte equality for both real core sparse rules\n'
    (out/(version+'.log')).write_text(text,encoding='utf-8',newline='\n')
    rows.append(dict(version=version,source_sha256=hashlib.sha256(source).hexdigest(),tests_run=1,failures=len(result.failures),errors=len(result.errors),skipped=len(result.skipped),windows_relative_paths_observed=len(seen)))
(out/'result.json').write_text(json.dumps(dict(host='Linux',windows_relative_path_emulation_only=True,actual_test_body_executed=True,old_failure_expected_and_retained=True,not_added_to_full_suite_counts=True,rows=rows),indent=2)+'\n')
print(rows)
