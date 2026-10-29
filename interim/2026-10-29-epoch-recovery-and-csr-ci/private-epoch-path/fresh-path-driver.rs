// Verification executable, not the production CLI. It links the exact private
// library after the focused suite passes and never rewrites solver identities.
use soil_fem_3d::{Error, Result, ground::Ground, mesh::Mesh,
    stages::{self, McEpochCheckpoint, McParameterMechanics, StageCase}};
use soil_fem_3d::mc_epoch_path::{Advance,Controller as PathController,Checkpoint as PathCheckpoint,Specification};
use std::{fs::{File,OpenOptions},io::{BufReader,Write}};

fn input<T:serde::de::DeserializeOwned>(path:&str)->Result<T> {
    Ok(serde_json::from_reader(BufReader::new(File::open(path)?))?)
}
fn output(value:&impl serde::Serialize,path:&str)->Result<()> {
    let bytes=serde_json::to_vec(value)?;
    let mut f=OpenOptions::new().write(true).create_new(true).open(path)?;
    f.write_all(&bytes)?;f.sync_all()?;Ok(())
}
fn c(controller:&McParameterMechanics<'_>,value:&str)->Result<Vec<soil_fem_3d::mc::Parameters>> {
    let value=value.parse::<f64>().map_err(|_|Error::Input("invalid cohesion".into()))?;
    Ok(controller.parameters().iter().cloned().map(|mut p| {p.cohesion_kpa=value;p}).collect())
}
fn run()->Result<()> {
    let a=std::env::args().collect::<Vec<_>>();
    if a.len()<3 {return Err(Error::Input("verification driver arguments missing".into()));}
    let mesh:Mesh=input(&a[2])?;let mesh=mesh.elevate_straight()?;
    match (a[1].as_str(),a.len()) {
        ("start-path",7)=>{
            let origin:stages::Run=input(&a[3])?;let spec:Specification=input(&a[4])?;
            let n=a[5].parse::<usize>().map_err(|_|Error::Input("invalid attempt count".into()))?;
            if n>10000 {return Err(Error::Input("attempt count exceeds verification cap".into()));}
            let mut path=PathController::new(&mesh,&origin,spec)?;
            for _ in 0..n {if path.advance()?==Advance::Finished {break;}}
            output(&path.checkpoint()?,&a[6])
        },
        ("resume-path",7)=>{
            let origin:stages::Run=input(&a[3])?;let cp:PathCheckpoint=input(&a[4])?;
            let n=a[5].parse::<usize>().map_err(|_|Error::Input("invalid attempt count".into()))?;
            if n>10000 {return Err(Error::Input("attempt count exceeds verification cap".into()));}
            let mut path=PathController::resume(&mesh,&origin,&cp)?;
            for _ in 0..n {if path.advance()?==Advance::Finished {break;}}
            output(&path.checkpoint()?,&a[6])
        },
        ("origin",6)=>{
            let ground:Ground=input(&a[3])?;let case:StageCase=input(&a[4])?;
            let run=stages::solve(&mesh,ground,case)?;
            if !run.complete {return Err(Error::Solve(run.failure_reason.unwrap_or_default()));}
            output(&run,&a[5])
        },
        ("advance",6)=>{
            let origin:stages::Run=input(&a[3])?;let mut m=McParameterMechanics::new(&mesh,&origin)?;
            m.try_parameters(c(&m,&a[4])?)?;
            output(&m.checkpoint()?,&a[5])
        },
        ("resume",7)=>{
            let origin:stages::Run=input(&a[3])?;let cp:McEpochCheckpoint=input(&a[4])?;
            let mut m=McParameterMechanics::resume(&mesh,&origin,&cp)?;
            let step=m.try_parameters(c(&m,&a[5])?)?;output(&step,&a[6])
        },
        ("continuous",7)=>{
            let origin:stages::Run=input(&a[3])?;let mut m=McParameterMechanics::new(&mesh,&origin)?;
            m.try_parameters(c(&m,&a[4])?)?;
            let step=m.try_parameters(c(&m,&a[5])?)?;output(&step,&a[6])
        },
        ("restore",6)=>{
            let origin:stages::Run=input(&a[3])?;let cp:McEpochCheckpoint=input(&a[4])?;
            let m=McParameterMechanics::resume(&mesh,&origin,&cp)?;output(&m.checkpoint()?,&a[5])
        },
        _=>Err(Error::Input("unknown verification operation/argument count".into())),
    }
}
fn main() {
    // This is deliberately a plain fresh process using the same exact image
    // for origin generation, continuation and recovery.
    if let Err(e)=run() {eprintln!("{e}");std::process::exit(1);}
}
