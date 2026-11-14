"""Build the existing read-only MCC observer in a NEW isolated crate."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]


def run(output):
    output.mkdir(parents=True,exist_ok=False);crate=output/'crate';crate.mkdir()
    shutil.copytree(ROOT/'src',crate/'src')
    for name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):shutil.copyfile(ROOT/name,crate/name)
    source={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'src').rglob('*.rs')}
    prior=ROOT/'interim/2026-11-09-nc-turning-path-diagnosis'
    for name,addition in [('mcc.rs','mcc-diagnostic-addition.rs'),('mcc_midpoint.rs','midpoint-diagnostic-addition.rs')]:
        p=crate/'src'/name;p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n')+(prior/addition).read_bytes())
    (crate/'examples').mkdir();shutil.copyfile(prior/'material_diagnostic.rs',crate/'examples/material_diagnostic.rs')
    subprocess.run(['cargo','fmt','--manifest-path',str(crate/'Cargo.toml')],check=True)
    subprocess.run(['cargo','build','--locked','--release','--manifest-path',str(crate/'Cargo.toml'),
                    '--target-dir',str(output/'target'),'--example','material_diagnostic'],check=True)
    binary=output/'target/release/examples'/('material_diagnostic.exe' if sys.platform=='win32' else 'material_diagnostic')
    if any(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=sha for p,sha in source.items()):raise ValueError('production source changed during observer preparation')
    record=dict(root_numerical_source_sha256=source,observer_binary=str(binary),observer_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                observer_only=True,production_modified=False)
    (output/'observer.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(run(p.parse_args().output.resolve()),indent=2))
