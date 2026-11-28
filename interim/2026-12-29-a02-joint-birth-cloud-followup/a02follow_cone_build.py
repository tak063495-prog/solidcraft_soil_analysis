"""Isolated projection observation, release vs dev builds; no production mutation."""
import argparse,difflib,json,os,shutil,subprocess
from pathlib import Path
from a02follow_run import ROOT,HERE,write,sha

def build(out):
    out.mkdir(parents=True,exist_ok=False);src=out/'source';src.mkdir()
    shutil.copytree(ROOT/'src',src/'src')
    for n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):shutil.copy2(ROOT/n,src/n)
    path=src/'src/joints_birth.rs';old=path.read_text();needle='        for di in selected {'
    assert old.count(needle)==1
    new=old.replace(needle,'''        if let Some(path) = std::env::var_os("A02FOLLOW_RAW_PROJECTION") {
            std::fs::write(path, serde_json::to_vec(&traction)?)?;
        }
'''+needle)
    path.write_text(new,encoding='utf-8',newline='\n')
    (out/'projection-observation.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True))))
    (src/'src/bin').mkdir(exist_ok=True);shutil.copy2(HERE/'a02follow_cone_probe.rs',src/'src/bin/a02follow_cone_probe.rs')
    manifest=(ROOT/'Cargo.toml').read_text().split('[lints.rust]')[0]
    (src/'Cargo.toml').write_text(manifest+'\n[lints.rust]\nunsafe_code = "forbid"\n')
    records=[]
    for profile,args in [('release',['--release']),('debug',[])]:
        cmd=['cargo','build','--locked',*args,'--bin','a02follow_cone_probe']
        with (out/(profile+'-build.log')).open('x') as log:r=subprocess.run(cmd,cwd=src,stdout=log,stderr=subprocess.STDOUT)
        exe=src/'target'/profile/('a02follow_cone_probe'+('.exe' if os.name=='nt' else ''))
        records.append(dict(profile=profile,command=cmd,exit_code=r.returncode,binary_sha256=sha(exe) if exe.exists() else None));r.check_returncode()
    write(out/'binding.json',dict(production_source_changed=False,builds=records,source_sha256={p.relative_to(src).as_posix():sha(p) for p in sorted((src/'src').rglob('*')) if p.is_file()}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);build(p.parse_args().output.resolve())
