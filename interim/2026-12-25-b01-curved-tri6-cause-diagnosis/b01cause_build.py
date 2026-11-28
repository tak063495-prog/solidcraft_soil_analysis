"""Build an isolated counterfactual from SHA-bound current element sources."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def generate(output):
    output.mkdir(parents=True, exist_ok=True)
    paths = ['src/plates.rs', 'src/plates/mindlin.rs', 'src/plates/offset.rs']
    sources = {p: (ROOT/p).read_bytes() for p in paths}
    original = sources[paths[1]].decode()
    start = original.index('        for (e, (i, j)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {')
    end = original.index('        let xy = x.map', start)
    removed = original[start:end]
    if '64. * f64::EPSILON * scale' not in removed or removed.count('return Err') != 1:
        raise ValueError('current geometry guard changed; refuse unreviewed transplant')
    # scale was used only by the removed guard. Keep it evaluated and mark it used.
    internal = original[:start]+'        let _ = scale; // original guard measurement retained\n'+original[end:]
    outer = sources[paths[0]].decode()
    outer = outer.replace('pub mod mindlin;', 'pub mod mindlin {\n'+internal+'\n}', 1)
    outer = outer.replace('pub mod offset;', 'pub mod offset {\n'+sources[paths[2]].decode()+'\n}', 1)
    generated = '#[allow(dead_code)]\nmod b01cause_flat {\n'+outer+'\n}\n'
    (output/'b01cause_flat.rs').write_text(generated, encoding='utf-8', newline='\n')
    shutil.copyfile(HERE/'b01cause_probe.rs', output/'b01cause_probe.rs')
    (output/'guard-only.diff').write_text(''.join(difflib.unified_diff(
        original.splitlines(True), internal.splitlines(True), fromfile=paths[1],
        tofile='work-only/b01cause_internal_mindlin.rs')), encoding='utf-8')
    for path, raw in sources.items():
        dest = output/'original'/path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(raw)
    binding = dict(diagnostic_only=True, production_guard_unchanged=True,
                   source_sha256={p: sha(b) for p, b in sources.items()},
                   generated_sha256=sha(generated.encode()), removed_guard_sha256=sha(removed.encode()),
                   modification='one straight-midpoint guard block omitted in work-only module; no physics changes')
    (output/'transplant-binding.json').write_text(json.dumps(binding, indent=2)+'\n')
    return binding


def build(output):
    binding = generate(output)
    deps = ROOT/'target/release/deps'
    def lib(name):
        files = list(deps.glob('lib'+name+'-*.rlib'))
        if not files:
            raise ValueError('run cargo build --locked --release first')
        return max(files, key=lambda p: p.stat().st_mtime)
    exe = output/('b01cause_probe.exe' if __import__('os').name == 'nt' else 'b01cause_probe')
    command = ['rustc', '--edition=2024', '-O', '-D', 'warnings',
               '--extern', 'soil_fem_3d='+str(ROOT/'target/release/libsoil_fem_3d.rlib'),
               '--extern', 'nalgebra='+str(lib('nalgebra')),
               '--extern', 'serde='+str(lib('serde')),
               '--extern', 'serde_json='+str(lib('serde_json')),
               '-L', 'dependency='+str(deps), str(output/'b01cause_probe.rs'), '-o', str(exe)]
    with (output/'build.log').open('w', encoding='utf-8') as log:
        p = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    binding.update(build_command=command, build_exit_code=p.returncode)
    if not p.returncode:
        binding['binary_sha256'] = sha(exe.read_bytes())
    (output/'build-binding.json').write_text(json.dumps(binding, indent=2)+'\n')
    p.check_returncode()
    return exe


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    print(build(args.output.resolve()))
