"""Link a read-only local probe to the unmodified public production kernel."""
import argparse, hashlib, json, os
from pathlib import Path
import subprocess
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);out=p.parse_args().output.resolve();out.mkdir(parents=True,exist_ok=False)
    lib=ROOT/'target/release/libsoil_fem_3d.rlib';deps=ROOT/'target/release/deps'
    serde=max(deps.glob('libserde_json-*.rlib'),key=lambda p:p.stat().st_mtime)
    exe=out/('a01follow_contact_probe.exe' if os.name=='nt' else 'a01follow_contact_probe')
    args=['rustc','--edition=2024','-O','-D','warnings','--extern','soil_fem_3d='+str(lib),'--extern','serde_json='+str(serde),
        '-L','dependency='+str(deps),str(HERE/'a01follow_contact_probe.rs'),'-o',str(exe)]
    with (out/'build.log').open('x') as log:result=subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    record=dict(command=args,exit_code=result.returncode,production_source_edited=False,
        sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (lib,HERE/'a01follow_contact_probe.rs',ROOT/'src/joint_law.rs')})
    if not result.returncode:record['binary_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
    (out/'binding.json').write_text(json.dumps(record,indent=2)+'\n');result.check_returncode()

if __name__=='__main__':main()
