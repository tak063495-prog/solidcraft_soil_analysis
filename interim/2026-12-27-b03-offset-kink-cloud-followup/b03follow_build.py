"""Link an explicit probe against unmodified production public kernels."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def build(output):
    output.mkdir(parents=True, exist_ok=True)
    deps = ROOT/'target/release/deps'
    args = ['rustc', '--edition=2024', '-O', '-D', 'warnings',
            '--extern', 'soil_fem_3d='+str(ROOT/'target/release/libsoil_fem_3d.rlib')]
    for name in ('nalgebra', 'serde_json'):
        lib = max(deps.glob('lib'+name+'-*.rlib'), key=lambda p: p.stat().st_mtime)
        args += ['--extern', name+'='+str(lib)]
    exe = output/('b03follow_probe.exe' if __import__('os').name == 'nt' else 'b03follow_probe')
    args += ['-L', 'dependency='+str(deps), str(HERE/'b03follow_probe.rs'), '-o', str(exe)]
    with (output/'build.log').open('w', encoding='utf-8') as log:
        p = subprocess.run(args, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    paths = ['src/plates.rs', 'src/plates/mindlin.rs', 'src/plates/offset.rs', 'src/beams.rs', 'tools/plate6_offset_reference.py', 'tools/plate6_reference.py', 'tools/plate_reference.py']
    binding = dict(command=args, exit_code=p.returncode, no_transplanted_source=True,
                   source_sha256={v: hashlib.sha256((ROOT/v).read_bytes()).hexdigest() for v in paths})
    if not p.returncode:
        binding['binary_sha256'] = hashlib.sha256(exe.read_bytes()).hexdigest()
    (output/'build-binding.json').write_text(json.dumps(binding, indent=2)+'\n')
    p.check_returncode()
    return exe


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    print(build(p.parse_args().output.resolve()))
