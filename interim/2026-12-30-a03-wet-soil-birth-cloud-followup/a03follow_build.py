"""Explicit isolated experiment; production sources are only read/copied."""
import argparse
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
from a03follow_binding import ROOT


def replace(text, old, new):
    if text.count(old) != 1:
        raise ValueError(('source replacement count', text.count(old), old))
    return text.replace(old, new)


def build(out):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = out/'source'
    source.mkdir()
    shutil.copytree(ROOT/'src', source/'src')
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
        shutil.copy2(ROOT/name, source/name)
    before = (source/'src/stages.rs').read_text(encoding='utf-8')
    text = replace(before, '        if dofs.transform.is_some() {', '''        for count in &mut partitions[self.point_count..] {
            *count = (*count).max(a03follow_minimum());
        }
        if dofs.transform.is_some() {''')
    text += '''
// Only the explicit isolated diagnostic build reads this switch.
fn a03follow_minimum() -> usize {
    std::env::var("A03FOLLOW_MINIMUM").ok().map(|v| v.parse().unwrap()).unwrap_or(1)
}
'''
    (source/'src/stages.rs').write_text(text, encoding='utf-8', newline='\n')
    (out/'stages.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True), text.splitlines(True))), encoding='utf-8')
    before = (source/'src/elements.rs').read_text(encoding='utf-8')
    text = replace(before, '    Ok(result)\n}\n\n#[cfg(test)]', '''    if std::env::var_os("A03FOLLOW_REVERSE_FACE_POINTS").is_some() {
        result.reverse();
    }
    Ok(result)
}

#[cfg(test)]''')
    (source/'src/elements.rs').write_text(text, encoding='utf-8', newline='\n')
    (out/'elements.patch').write_text(''.join(difflib.unified_diff(before.splitlines(True), text.splitlines(True))), encoding='utf-8')
    (source/'src/bin').mkdir(exist_ok=True)
    shutil.copy2(ROOT/'examples/wet_birth_a03_probe.rs', source/'src/bin/a03follow_probe.rs')
    env = dict(os.environ, CARGO_TARGET_DIR=str(out/'target'))
    commands = []
    for args in (['cargo', 'fmt', '--all'], ['cargo', 'build', '--locked', '--release', '--bin', 'a03follow_probe']):
        with (out/('fmt.log' if 'fmt' in args else 'build.log')).open('x') as log:
            p = subprocess.run(args, cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(command=args, exit_code=p.returncode))
        if p.returncode:
            raise RuntimeError(commands)
    exe = out/('target/release/a03follow_probe'+('.exe' if os.name == 'nt' else ''))
    value = dict(base='89e9dfc26cbe444199179d7faa1edb1437849309', commands=commands,
                 production_edited=False, executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                 source_sha256={p.relative_to(source).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(source.rglob('*')) if p.is_file()})
    (out/'binding.json').write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    print(exe)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    build(p.parse_args().output)
