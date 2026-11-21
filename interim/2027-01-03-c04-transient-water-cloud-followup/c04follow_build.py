"""Isolated P/R exchange inherited exactly from A-04; observational metadata only."""
import argparse
import difflib
import os
import shutil
import subprocess
from c04follow_common import ROOT, BASE, load, sha, write

A04 = ROOT/'interim/2026-12-31-a04-richards-wet-joint-cloud-followup'

def stage_patch(text):
    import sys
    sys.path.insert(0, str(A04))
    prior = load('c04follow_inherited_build', A04/'a04follow_build.py')
    return prior.patch_stages(text).replace('a04follow', 'c04follow').replace('A04FOLLOW', 'C04FOLLOW')

def replace(text, before, after):
    if text.count(before) != 1:
        raise ValueError(('source anchor mismatch', text.count(before), before[:100]))
    return text.replace(before, after)

def law_patch(text):
    text = replace(text, 'pub struct PointState {', '''pub struct PointState {
    // Isolated observation only. Never read by constitutive or solver decisions.
    #[serde(default)]
    pub c04follow_decision: Option<serde_json::Value>,''')
    text = replace(text, '        Ok(Self {\n            birth: None,', '''        Ok(Self {
            c04follow_decision: None,
            birth: None,''')
    text = replace(text, '        let q = normal * n + tau;', '''        s.c04follow_decision = Some(serde_json::json!({
            "kind": "closed", "trial_norm_kpa": norm, "radius_kpa": limit,
            "flow_m": <[f64; 3]>::from(plastic), "gap_m": (target-old.elastic_offset()).dot(&n),
            "status_lhs": plastic.norm(), "status_rhs": 0., "status_operator": ">",
            "return_shear_norm_kpa": tau.norm(), "returned_margin_kpa": limit-tau.norm()
        }));
        let q = normal * n + tau;''')
    # On any open endpoint this metadata is explicitly marked, never inherited as a closed decision.
    text = text.replace('            s.status = Status::Open;', '''            s.c04follow_decision = Some(serde_json::json!({"kind":"open", "gap_m":g1}));
            s.status = Status::Open;''')
    text = replace(text, '            return Ok(old.clone());', '''            let mut observed = old.clone();
            observed.c04follow_decision = Some(serde_json::json!({"kind":"unchanged", "inherited_status":old.status}));
            return Ok(observed);''')
    return text

def build(output):
    output.mkdir(parents=True, exist_ok=False)
    source = output/'source'; source.mkdir()
    shutil.copytree(ROOT/'src', source/'src')
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
        shutil.copy2(ROOT/name, source/name)
    for name, transform in (('stages.rs', stage_patch), ('joint_law.rs', law_patch)):
        path = source/'src'/name
        before = path.read_text(encoding='utf-8'); after = transform(before)
        path.write_text(after, encoding='utf-8', newline='\n')
        (output/(name+'.patch')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True))),encoding='utf-8')
    binary = source/'src/bin/c04follow_probe.rs'; binary.parent.mkdir(exist_ok=True)
    probe = (ROOT/'examples/wet_joint_transient_probe.rs').read_text()
    probe = replace(probe, '    assert!(mesh.tetrahedra.len() <= 768 || std::env::var_os("C04_USER_PC").is_some());',
                    '    assert_eq!(mesh.tetrahedra.len(), 12);\n    assert_eq!(mesh.nodes.len(), 16);')
    binary.write_text(probe, encoding='utf-8', newline='\n')
    env = dict(os.environ, CARGO_TARGET_DIR=str(output/'target'))
    commands = []
    for name,args in (('fmt',['cargo','fmt','--all']),('build',['cargo','build','--locked','--release','--bin','c04follow_probe'])):
        with (output/(name+'.log')).open('x') as log:
            p = subprocess.run(args, cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(command=args, exit_code=p.returncode))
        if p.returncode:
            write(output/'build-failure.json',dict(complete=False, commands=commands)); raise RuntimeError(commands)
    suffix = '.exe' if os.name == 'nt' else ''
    write(output/'binding.json',dict(base=BASE,production_edited=False,commands=commands,
        probe_sha256=sha(output/('target/release/c04follow_probe'+suffix)),
        sources={p.relative_to(source).as_posix():sha(p) for p in sorted(source.rglob('*')) if p.is_file()}))

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=__import__('pathlib').Path)
    build(p.parse_args().output.resolve())
