"""Build a public-API observer with the unchanged locked production dependency."""
import argparse
import json
import subprocess
from pathlib import Path
from common import ROOT, HERE

def build(destination):
    destination.mkdir(parents=True, exist_ok=False)
    (destination/'Cargo.toml').write_text('[package]\nname="b02-residual-observer"\nversion="0.0.0"\nedition="2024"\n[dependencies]\nsoil-fem-3d={path='+json.dumps(str(ROOT))+'}\nserde={version="1",features=["derive"]}\nserde_json={version="1",features=["float_roundtrip"]}\n[[bin]]\nname="b02-residual-observer"\npath='+json.dumps(str(HERE/'observer.rs'))+'\n')
    (destination/'Cargo.lock').write_text((ROOT/'Cargo.lock').read_text()+'\n[[package]]\nname="b02-residual-observer"\nversion="0.0.0"\ndependencies=["serde","serde_json","soil-fem-3d"]\n')
    for args in (['fmt','--manifest-path',str(destination/'Cargo.toml'),'--check'],
                 ['clippy','--locked','--release','--manifest-path',str(destination/'Cargo.toml'),'--target-dir',str(ROOT/'target'),'--','-D','warnings'],
                 ['build','--locked','--release','--manifest-path',str(destination/'Cargo.toml'),'--target-dir',str(ROOT/'target')]):
        subprocess.run(['cargo',*args], cwd=ROOT, check=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination',type=Path);a=p.parse_args();build(a.destination.resolve())
