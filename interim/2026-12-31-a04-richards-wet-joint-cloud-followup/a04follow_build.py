"""Build isolated opt-in controls; never edits production Rust files."""
import argparse
import difflib
import os
from pathlib import Path
import shutil
import subprocess
from a04follow_common import ROOT, BASE, write, sha

def replace(t,a,b):
    if t.count(a)!=1:
        raise ValueError(('source anchor mismatch',t.count(a),a[:100]))
    return t.replace(a,b)

def patch_stages(t):
    t=replace(t,'        let diagnostic = self\n', '''        // This switch exists only in the isolated follow-up source copy.
        let a04follow_p = std::env::var("A04FOLLOW_P").unwrap_or("N".into());
        let a04follow_r = std::env::var("A04FOLLOW_R").unwrap_or("N".into());
        assert!(["N", "I"].contains(&a04follow_p.as_str()));
        assert!(["N", "I"].contains(&a04follow_r.as_str()));
        if a04follow_p == "I" {
            partitions.fill(1);
        }
        let diagnostic = self
''')
    t=replace(t,'            if relative <= settings.relative_force_tolerance {\n', '''            if relative <= settings.relative_force_tolerance {
                if a04follow_r == "I" {
                    let needed = self.predictor_partitions(old, &u, &active, settings)?;
                    let mut changed = false;
                    for (count, n) in partitions.iter_mut().zip(needed) {
                        if n > *count {
                            *count = n;
                            changed = true;
                        }
                    }
                    if changed { continue; }
                }
''')
    # The small input has supported bodies even when the interface is open.
    # Fail rather than silently claiming exchange of a neutral/cold predictor.
    t=replace(t,'            let cold = if iteration == 0 && gauges.is_some() {', '''            if (a04follow_p == "I" || a04follow_r == "I") && gauges.is_some() {
                return Err(Error::Control("a04follow exchange unavailable for neutral/cold predictors".into()));
            }
            let cold = if iteration == 0 && gauges.is_some() {''')
    t=replace(t,'            for _ in 0..16 {\n                let mut candidate = u.clone();', '''            for _ in 0..if a04follow_p == "I" { 18 } else { 16 } {
                let mut candidate = u.clone();''')
    a='                for (ei, e) in self.elements.iter().enumerate() {\n                    if active[ei] {\n                        let du = DVector::from_fn'
    b='                if a04follow_r == "N" {\n'+a
    split=t.index('            let mut alpha = 1.;',t.index('    fn trial_with_predictor('))
    t=t[:split]+replace(t[split:],a,b)
    t=replace(t,'                if repartition {\n                    if let Some(o) = diagnostic {','                }\n                if repartition {\n                    if let Some(o) = diagnostic {')
    t=replace(t,'                    if self.reduce(&r, &dofs).norm() < reduced.norm() * (1. - 1e-4 * alpha) {', '''                    if self.reduce(&r, &dofs).norm() < reduced.norm() * (1. - if a04follow_p == "I" {0.} else {1e-4} * alpha) {''')
    return t

def build(out):
    out.mkdir(parents=True,exist_ok=False); source=out/'source';source.mkdir()
    shutil.copytree(ROOT/'src',source/'src')
    for name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml'):
        shutil.copy2(ROOT/name,source/name)
    for filename,transform in (
        ('stages.rs',patch_stages),
        ('elements.rs',lambda t:replace(t,'    Ok(result)\n}\n\n#[cfg(test)]', '''    if std::env::var_os("A04FOLLOW_Q").is_some() { result.reverse(); }
    Ok(result)
}

#[cfg(test)]''')),
        ('mixed_flow.rs',lambda t:replace(t,'        let volume_m3: f64 = points.iter().map(|p| p.weight_m3).sum();', '''        if std::env::var_os("A04FOLLOW_Q").is_some() { points.reverse(); }
        let volume_m3: f64 = points.iter().map(|p| p.weight_m3).sum();'''))):
        path=source/'src'/filename;before=path.read_text(encoding='utf-8');after=transform(before)
        path.write_text(after,encoding='utf-8',newline='\n')
        (out/(filename+'.patch')).write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True))),encoding='utf-8')
    (source/'src/bin').mkdir(exist_ok=True)
    for target,example in (('a04follow_c04_probe','wet_joint_transient_probe'),('a04follow_dry_probe','joint_a02_history_probe'),('a04follow_hydro_probe','wet_birth_a03_probe')):
        shutil.copy2(ROOT/('examples/'+example+'.rs'),source/('src/bin/'+target+'.rs'))
    env=dict(os.environ,CARGO_TARGET_DIR=str(out/'target'))
    commands=[]
    for name,args in (('fmt',['cargo','fmt','--all']),('build',['cargo','build','--locked','--release','--bins'])):
        with (out/(name+'.log')).open('x') as f:
            p=subprocess.run(args,cwd=source,env=env,stdout=f,stderr=subprocess.STDOUT)
        commands.append(dict(command=args,exit_code=p.returncode))
        if p.returncode:raise RuntimeError(commands)
    suffix='.exe' if os.name=='nt' else ''
    binaries={n:sha(out/('target/release/'+n+suffix)) for n in ('soil-fem-3d','a04follow_c04_probe','a04follow_dry_probe','a04follow_hydro_probe')}
    write(out/'binding.json',dict(base=BASE,production_edited=False,commands=commands,binaries=binaries,
          sources={p.relative_to(source).as_posix():sha(p) for p in sorted(source.rglob('*')) if p.is_file()}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);build(p.parse_args().output.resolve())
