"""Build the public-API diagnostic in work/, outside production Cargo targets."""
import argparse
import json
import os
from pathlib import Path
import subprocess
from cases import ROOT, write


def build(destination):
    destination.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).with_name('observer.rs').resolve()
    manifest = destination/'Cargo.toml'
    manifest.write_text('[package]\nname = "b02-observer"\nversion = "0.0.0"\nedition = "2024"\n'
                        '[dependencies]\nsoil-fem-3d = { path = '+json.dumps(str(ROOT))+' }\n'
                        'serde = { version = "1", features = ["derive"] }\nserde_json = { version = "1", features = ["float_roundtrip"] }\n'
                        '[[bin]]\nname = "b02-observer"\npath = '+json.dumps(str(source))+'\n', encoding='utf-8')
    # Preserve every pinned production dependency. Only the diagnostic root
    # package is appended; cargo --locked verifies the resulting graph.
    (destination/'Cargo.lock').write_text((ROOT/'Cargo.lock').read_text()+
        '\n[[package]]\nname = "b02-observer"\nversion = "0.0.0"\ndependencies = ["serde", "serde_json", "soil-fem-3d"]\n', encoding='utf-8')
    subprocess.run(['cargo', 'build', '--locked', '--release', '--manifest-path', str(manifest), '--target-dir', str(ROOT/'target')], check=True, cwd=ROOT)
    return ROOT/'target/release'/('b02-observer.exe' if os.name == 'nt' else 'b02-observer')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('destination', type=Path); a=p.parse_args()
    print(build(a.destination.resolve()))
