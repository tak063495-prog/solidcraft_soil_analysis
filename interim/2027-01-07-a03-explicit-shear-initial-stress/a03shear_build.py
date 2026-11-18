"""Build an unmodified production solver copy and the existing observation probe."""
import argparse, hashlib, json, os
from pathlib import Path
import shutil, subprocess
ROOT=Path(__file__).resolve().parents[2]

def build(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);source=out/'source';source.mkdir()
    shutil.copytree(ROOT/'src',source/'src')
    # Include compile-time test fixtures, without introducing more example targets.
    (source/'examples').mkdir()
    for p in (ROOT/'examples').glob('*.json'):shutil.copy2(p,source/'examples'/p.name)
    for n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):shutil.copy2(ROOT/n,source/n)
    (source/'src/bin').mkdir(exist_ok=True)
    shutil.copy2(ROOT/'examples/wet_birth_a03_probe.rs',source/'src/bin/a03shear_probe.rs')
    (out/'stages.patch').write_text('');(out/'elements.patch').write_text('')
    env=dict(os.environ,CARGO_TARGET_DIR=str(out/'target'),CARGO_BUILD_JOBS='2')
    commands=[]
    for name,args in [('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('build',['cargo','build','--locked','--release','--bin','a03shear_probe'])]:
        with (out/(name+'.log')).open('x') as log:p=subprocess.run(args,cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT)
        commands.append(dict(command=args,exit_code=p.returncode));p.check_returncode()
    exe=out/('target/release/a03shear_probe'+('.exe' if os.name=='nt' else ''))
    (out/'binding.json').write_text(json.dumps(dict(base='cdaf850c0347c8d841decce7ddff608d03ae7933',commands=commands,production_edited=False,
        copied_solver_sources_exact=all((source/p.relative_to(ROOT)).read_bytes()==p.read_bytes() for p in (ROOT/'src').rglob('*.rs')),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),source_sha256={p.relative_to(source).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.rglob('*')) if p.is_file()}),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);build(p.parse_args().output)
