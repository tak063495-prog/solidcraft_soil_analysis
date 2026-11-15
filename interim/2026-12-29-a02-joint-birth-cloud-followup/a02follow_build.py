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
    text=replace(text,'        if self.joints.is_none() {','        if a02follow_policy("P") == "independent" || self.joints.is_none() {')
    text=replace(text,'        let contact_trust = if let Some(geometry)', '        let contact_trust = if a02follow_policy("P") == "independent" { false } else if let Some(geometry)')
    text=replace(text,'            partitions[self.point_count..].copy_from_slice(&needed);','            if a02follow_policy("R") != "independent" { partitions[self.point_count..].copy_from_slice(&needed); }')
    text=replace(text,'            if relative <= settings.relative_force_tolerance {','''            if relative <= settings.relative_force_tolerance {
                if a02follow_policy("R") == "independent" {
                    let required = self.predictor_partitions(old, &u, &active, settings)?;
                    let mut changed = false;
                    for i in self.point_count..partitions.len() {
                        if required[i] > partitions[i] {partitions[i] = required[i];changed = true;}
                    }
                    if changed {continue;}
                }''')
    text=replace(text,'            let gauges = if let Some(geometry)', '            let gauges = if a02follow_policy("P") == "independent" { None } else if let Some(geometry)')
    start=text.index('                for (ei, e) in self.elements.iter().enumerate() {',text.index('            for _ in 0..16 {'))
    end=text.index('                if repartition {',start)
    text=text[:start]+'                if a02follow_policy("R") != "independent" {\n'+text[start:end]+'                }\n'+text[end:]
    text+='''
// No production caller or environment variable reads these diagnostic policies.
fn a02follow_policy(name: &str) -> String {
    std::env::var(format!("A02FOLLOW_{name}")).unwrap_or_else(|_| "native".into())
}
'''
    text=replace(text,'        if dofs.transform.is_some() {','''        for count in &mut partitions[self.point_count..] { *count = (*count).max(a02follow_minimum()); }
        if dofs.transform.is_some() {''')
    text += '\nfn a02follow_minimum() -> usize { std::env::var("A02FOLLOW_MINIMUM").ok().map(|v| v.parse().unwrap()).unwrap_or(1) }\n'
    text=replace(text,'            let relative = reduced.norm() / scale;','''            let relative = reduced.norm() / scale;
            if let Some(o) = self.observer.filter(|_| std::env::var_os("A02FOLLOW_TRACE").is_some()) {
                o.on_diagnostic(&serde_json::json!({"kind":"a02follow_iterate", "stage":stage, "progress":progress,
                    "iteration":iteration,"partitions":partitions,"relative_force_residual":relative,"projected_load":projected_load}))?;
            }''')
    text=replace(text,'                if repartition {','''                if repartition {
                    if let Some(o) = self.observer.filter(|_| std::env::var_os("A02FOLLOW_TRACE").is_some()) {
                        o.on_diagnostic(&serde_json::json!({"kind":"a02follow_repartition","stage":stage,"progress":progress,"iteration":iteration,"partitions":partitions,"alpha":alpha}))?;
                    }''')
    (source/'src/stages.rs').write_text(text)
    (out/'stages.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),text.splitlines(True),fromfile='production/src/stages.rs',tofile='isolated/src/stages.rs')))
    shutil.copy2(HERE/'a02follow_service.rs',source/'src/a02follow_service.rs')
    with (source/'src/lib.rs').open('a') as f:f.write('\npub mod a02follow_service;\n')
    probe=(ROOT/'examples/joint_a02_history_probe.rs').read_text()
    probe=replace(probe,'    diagnostic: bool,','    diagnostic: bool,\n    events: RefCell<BufWriter<std::fs::File>>,')
    probe=replace(probe,'        if value["kind"] == "iterate" {','''        let mut events = self.events.borrow_mut();
        let mut compact = serde_json::Map::new();
        for key in ["kind", "stage", "progress", "iteration", "partitions", "relative_force_residual", "alpha", "projected_load", "error"] {
            if let Some(v) = value.get(key) { compact.insert(key.into(), v.clone()); }
        }
        serde_json::to_writer(&mut *events, &compact)?; writeln!(events)?; events.flush()?;
        if false {''')
    probe=replace(probe,'        diagnostic: false,','''        diagnostic: false,
        events: RefCell::new(BufWriter::new(OpenOptions::new().write(true).create_new(true).open(folder.join("events.jsonl"))?)),''')
    (source/'src/bin').mkdir(exist_ok=True)
    (source/'src/bin/a02follow_probe.rs').write_text(probe)
    (source/'src/bin/a02follow_predictor.rs').write_text('''use std::io::{self, BufRead, Write};
fn main() {let stdin=io::stdin();let mut out=io::stdout().lock();for line in stdin.lock().lines() {
let value=serde_json::from_str(&line.unwrap()).unwrap();
let answer=match soil_fem_3d::a02follow_service::predictor(value) {Ok(v)=>v,Err(e)=>serde_json::json!({"error":e.to_string()})};
serde_json::to_writer(&mut out,&answer).unwrap();writeln!(out).unwrap();out.flush().unwrap();}}
''')
    # A standalone manifest avoids auto-discovering repository examples/tests.
    cargo=(ROOT/'Cargo.toml').read_text().split('[lints.rust]')[0]
    (source/'Cargo.toml').write_text(cargo+'\n[lints.rust]\nunsafe_code = "forbid"\n')
    target=Path(os.environ.get('A02FOLLOW_TARGET_DIR',str(out/'target'))).resolve()
    env=dict(os.environ,CARGO_TARGET_DIR=str(target))
    with (out/'build.log').open('w') as log:
        p=subprocess.run(['cargo','build','--locked','--release','--bin','a02follow_probe','--bin','a02follow_predictor'],cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT)
    binding=dict(exit_code=p.returncode,base='11df919ad053fc742fa40cb986953d8f53de28db',production_edited=False,
        source_sha256={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.rglob('*')) if p.is_file()})
    if not p.returncode:
        (out/'target/release').mkdir(parents=True,exist_ok=True)
        for name in ('a02follow_probe','a02follow_predictor'):
            exe=name+('.exe' if os.name=='nt' else '')
            if target!=out/'target':shutil.copy2(target/'release'/exe,out/'target/release'/exe)
        binding['binary_sha256']={name:hashlib.sha256((out/'target/release'/(name+('.exe' if os.name=='nt' else ''))).read_bytes()).hexdigest() for name in ('a02follow_probe','a02follow_predictor')}
    (out/'binding.json').write_text(json.dumps(binding,indent=2)+'\n');p.check_returncode()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);build(p.parse_args().output.resolve())
