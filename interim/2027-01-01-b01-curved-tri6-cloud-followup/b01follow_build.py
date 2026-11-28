"""Build only the observation caller; never copy or alter the element."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
from b01follow_common import BASE, HERE, ROOT, sha, write


def build(output):
    output.mkdir(parents=True, exist_ok=False)
    sources = ['src/plates.rs', 'src/plates/mindlin.rs', 'src/plates/offset.rs',
               'Cargo.toml', 'Cargo.lock', str(HERE.relative_to(ROOT)/'b01follow_probe.rs')]
    for path in sources:
        target = output/'source-snapshot'/path
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/path, target)
    binding = dict(base=BASE, diagnostic_only=True, production_guard_unchanged=True,
                   kernel_modifications=0, source_sha256={p:sha(ROOT/p) for p in sources})
    deps = ROOT/'target/release/deps'
    def lib(name):
        files = list(deps.glob('lib'+name+'-*.rlib'))
        if not files:
            raise ValueError('cargo build --locked --release is required')
        return max(files, key=lambda p:p.stat().st_mtime)
    exe = output/('b01follow_probe.exe' if os.name == 'nt' else 'b01follow_probe')
    command = ['rustc', '--edition=2024', '-O', '-D', 'warnings', '--extern',
               'soil_fem_3d='+str(ROOT/'target/release/libsoil_fem_3d.rlib'),
               '--extern', 'nalgebra='+str(lib('nalgebra')), '--extern', 'serde_json='+str(lib('serde_json')),
               '-L', 'dependency='+str(deps), str(HERE/'b01follow_probe.rs'), '-o', str(exe)]
    with (output/'build.log').open('x', encoding='utf-8') as log:
        process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    binding.update(command=command, exit_code=process.returncode)
    if not process.returncode:
        binding['binary_sha256'] = sha(exe)
    write(output/'build-binding.json', binding); process.check_returncode()
    return exe


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    print(build(p.parse_args().output.resolve()))
