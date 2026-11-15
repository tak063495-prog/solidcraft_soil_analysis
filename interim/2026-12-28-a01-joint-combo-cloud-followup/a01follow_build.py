"""Build an explicit, isolated policy experiment; never edit production source."""
import argparse, difflib, hashlib, json, os
from pathlib import Path
import shutil, subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def replace(text, old, new):
    if text.count(old)!=1:raise ValueError(('source replacement count',text.count(old),old[:100]))
    return text.replace(old,new)

def build(out):
    out.mkdir(parents=True,exist_ok=False)
    source=out/'source';source.mkdir()
    shutil.copytree(ROOT/'src',source/'src')
    for name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):shutil.copy2(ROOT/name,source/name)
    original=(source/'src/stages.rs').read_text()
    text=original
    text=replace(text,'        if self.joints.is_none() {','        if a01follow_policy("P") == "independent" || self.joints.is_none() {')
    text=replace(text,'        let contact_trust = if let Some(geometry)', '        let contact_trust = if a01follow_policy("P") == "independent" { false } else if let Some(geometry)')
    text=replace(text,'            partitions[self.point_count..].copy_from_slice(&needed);','            if a01follow_policy("R") != "independent" { partitions[self.point_count..].copy_from_slice(&needed); }')
    text=replace(text,'            if relative <= settings.relative_force_tolerance {','''            if relative <= settings.relative_force_tolerance {
                if a01follow_policy("R") == "independent" {
                    let required = self.predictor_partitions(old, &u, &active, settings)?;
                    let mut changed = false;
                    for i in self.point_count..partitions.len() {
                        if required[i] > partitions[i] {partitions[i] = required[i];changed = true;}
                    }
                    if changed {continue;}
                }''')
    text=replace(text,'            let gauges = if let Some(geometry)', '            let gauges = if a01follow_policy("P") == "independent" { None } else if let Some(geometry)')
    start=text.index('                for (ei, e) in self.elements.iter().enumerate() {',text.index('            for _ in 0..16 {'))
    end=text.index('                if repartition {',start)
    text=text[:start]+'                if a01follow_policy("R") != "independent" {\n'+text[start:end]+'                }\n'+text[end:]
    text+='''
// No production caller or environment variable reads these diagnostic policies.
fn a01follow_policy(name: &str) -> String {
    std::env::var(format!("A01FOLLOW_{name}")).unwrap_or_else(|_| "native".into())
}
'''
    (source/'src/stages.rs').write_text(text)
    (out/'stages.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='production/src/stages.rs',tofile='isolated/src/stages.rs')))
    shutil.copy2(HERE/'a01follow_service.rs',source/'src/a01follow_service.rs')
    with (source/'src/lib.rs').open('a') as f:f.write('\npub mod a01follow_service;\n')
    probe=(ROOT/'examples/joint_combo_probe.rs').read_text()
    probe=probe.replace('    assert!(run.complete);\n    println!', '    println!')
    # Preserve partial accepted prefix on numerical failure, instead of asserting complete.
    probe=replace(probe,'            assert!(run.complete);\n            break;', '            break;')
    probe=replace(probe,'        if value["kind"] == "iterate" {','''        {
            let mut stream = self.stream.borrow_mut();
            serde_json::to_writer(&mut *stream, value)?;
            writeln!(stream)?;stream.flush()?;
        }
        if value["kind"] == "iterate" {''')
    probe=probe.replace('    fn stop_requested(&self)', '    fn a01_joint_cause_enabled(&self) -> bool { true }\n    fn stop_requested(&self)')
    (source/'src/bin').mkdir(exist_ok=True)
    (source/'src/bin/a01follow_probe.rs').write_text(probe)
    (source/'src/bin/a01follow_predictor.rs').write_text('''use std::io::{self, BufRead, Write};
fn main() {let stdin=io::stdin();let mut out=io::stdout().lock();for line in stdin.lock().lines() {
let value=serde_json::from_str(&line.unwrap()).unwrap();
let answer=match soil_fem_3d::a01follow_service::predictor(value) {Ok(v)=>v,Err(e)=>serde_json::json!({"error":e.to_string()})};
serde_json::to_writer(&mut out,&answer).unwrap();writeln!(out).unwrap();out.flush().unwrap();}}
''')
    # A standalone manifest avoids auto-discovering repository examples/tests.
    cargo=(ROOT/'Cargo.toml').read_text().split('[lints.rust]')[0]
    (source/'Cargo.toml').write_text(cargo+'\n[lints.rust]\nunsafe_code = "forbid"\n')
    target=Path(os.environ.get('A01FOLLOW_TARGET_DIR',str(out/'target'))).resolve()
    env=dict(os.environ,CARGO_TARGET_DIR=str(target))
    with (out/'build.log').open('w') as log:
        p=subprocess.run(['cargo','build','--locked','--release','--bin','a01follow_probe','--bin','a01follow_predictor'],cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT)
    binding=dict(exit_code=p.returncode,base='b383aa007afd99e74f5988060335b2e07c8d341f',production_edited=False,
        source_sha256={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.rglob('*')) if p.is_file()})
    if not p.returncode:
        (out/'target/release').mkdir(parents=True,exist_ok=True)
        for name in ('a01follow_probe','a01follow_predictor'):
            exe=name+('.exe' if os.name=='nt' else '')
            if target!=out/'target':shutil.copy2(target/'release'/exe,out/'target/release'/exe)
        binding['binary_sha256']={name:hashlib.sha256((out/'target/release'/(name+('.exe' if os.name=='nt' else ''))).read_bytes()).hexdigest() for name in ('a01follow_probe','a01follow_predictor')}
    (out/'binding.json').write_text(json.dumps(binding,indent=2)+'\n');p.check_returncode()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);build(p.parse_args().output.resolve())
