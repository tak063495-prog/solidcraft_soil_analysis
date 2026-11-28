"""Reproduce Git path quoting directly, without PowerShell/cmd or FEM solves."""
import argparse, hashlib, json, subprocess, tempfile
from pathlib import Path, PureWindowsPath
ROOT=Path(__file__).resolve().parents[2]
ORIGINAL=ROOT/'interim/2027-01-06-a03-unequal-thickness-anisotropic-k'

def reproduce(out):
    out.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((ORIGINAL/'archive-index.json').read_text())
    paths=[(Path('interim')/ORIGINAL.name/p['name']).as_posix() for p in manifest['parts']]
    winpaths=[str(PureWindowsPath(p)) for p in paths]
    prior=json.loads((ROOT/'interim/2026-12-30-a03-wet-soil-birth-cloud-followup/validation/sparse-checkout-contract.json').read_text())
    rows=[]
    with tempfile.TemporaryDirectory() as temp:
        rules=Path(temp)/'rules';rules.write_bytes(('\n'.join(prior['actual_core_sparse_rules']['python-core-1']['rules'])+'\n').encode())
        for style,selected in (('posix',paths),('windows',winpaths)):
            for newline,separator in (('lf','\n'),('crlf','\r\n')):
                data=(separator.join(selected)+separator).encode()
                args=['git','sparse-checkout','check-rules','--no-cone','--rules-file',str(rules)]
                p=subprocess.run(args,cwd=ROOT,input=data,capture_output=True,check=True)
                name=style+'-'+newline
                (out/(name+'-stdin.bin')).write_bytes(data);(out/(name+'-stdout.bin')).write_bytes(p.stdout)
                rows.append(dict(name=name,command=args,shell=False,stdin_hex=data.hex(),stdout_hex=p.stdout.hex(),
                                 old_line_comparison_equal=p.stdout.decode().splitlines()==data.decode().splitlines()))
        data=b''.join(p.encode()+b'\0' for p in paths)
        p=subprocess.run(['git','sparse-checkout','check-rules','-z','--no-cone','--rules-file',str(rules)],cwd=ROOT,input=data,capture_output=True,check=True)
        assert p.stdout==data
        (out/'fixed-stdin.bin').write_bytes(data);(out/'fixed-stdout.bin').write_bytes(p.stdout)
    assert not next(r for r in rows if r['name']=='windows-lf')['old_line_comparison_equal']
    report=dict(passed=False,a03_complete=False,diagnostic_only=True,git_version=subprocess.check_output(['git','--version'],text=True).strip(),
                reproduction_only_not_windows_execution=True,old_windows_style_reproduces_mismatch_on_linux=True,
                fixed_posix_nul_bytes_exact=True,new_fem_solves=0,new_accepted_states=0,rows=rows,
                files={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.glob('*.bin'))})
    (out/'reproduction.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print({k:v for k,v in report.items() if k not in ('rows','files')})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);reproduce(p.parse_args().output)
