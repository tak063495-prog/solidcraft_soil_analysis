"""Build a separate diagnostic executable with the production locked graph."""
import argparse
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).with_name('observer.rs').resolve()
    (out/'Cargo.toml').write_text('[package]\nname="p602-observer"\nversion="0.0.0"\nedition="2024"\n'
        '[dependencies]\nsoil-fem-3d={path='+json.dumps(str(ROOT))+'}\n'
        'serde={version="1",features=["derive"]}\nserde_json={version="1",features=["float_roundtrip"]}\n'
        '[[bin]]\nname="p602-observer"\npath='+json.dumps(str(source))+'\n', encoding='utf-8')
    (out/'Cargo.lock').write_text((ROOT/'Cargo.lock').read_text()+
        '\n[[package]]\nname="p602-observer"\nversion="0.0.0"\ndependencies=["serde","serde_json","soil-fem-3d"]\n', encoding='utf-8')
    subprocess.run(['cargo','build','--locked','--release','--manifest-path',str(out/'Cargo.toml'),
                    '--target-dir',str(ROOT/'target')], cwd=ROOT, check=True)
    return ROOT/'target/release'/('p602-observer.exe' if os.name=='nt' else 'p602-observer')


if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('output',type=Path); a=p.parse_args()
    print(build(a.output.resolve()))
