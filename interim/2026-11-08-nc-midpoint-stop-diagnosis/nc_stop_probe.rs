//! Read-only, bounded replay of recorded trials, never a resumable Run.
use soil_fem_3d::{Error, Result, stages, stage_control::{Observer, Progress}};
use std::{cell::RefCell, fs::{self, File, OpenOptions}, io::Write, path::Path};
struct Trace(RefCell<File>);
impl Observer for Trace {
    fn on_progress(&self, p: &Progress) -> Result<()> {
        self.on_diagnostic(&serde_json::json!({"kind":"progress", "event":p}))
    }
    fn diagnostic_enabled(&self, _: usize, _: f64) -> bool { true }
    fn on_diagnostic(&self, v: &serde_json::Value) -> Result<()> {
        let mut f=self.0.borrow_mut();serde_json::to_writer(&mut *f,v)?;writeln!(f)?;Ok(())
    }
}
fn main() -> Result<()> {
    let a:Vec<_>=std::env::args().collect();
    if a.len()!=5 { return Err(Error::Input("nc_stop_probe MESH CHECKPOINT TARGET NEW_DIR".into())); }
    let mesh=serde_json::from_slice(&fs::read(&a[1])?)?;
    let saved=serde_json::from_slice(&fs::read(&a[2])?)?;
    let next=a[3].parse().map_err(|_|Error::Input("invalid target".into()))?;
    let out=Path::new(&a[4]);fs::create_dir(out)?;
    let trace=Trace(RefCell::new(OpenOptions::new().write(true).create_new(true).open(out.join("trace.jsonl"))?));
    let result=stages::diagnose_recorded_nc_trial(&mesh,&saved,next,&trace)?;
    fs::write(out.join("diagnostic.json"),serde_json::to_vec_pretty(&result)?)?;
    println!("{}",serde_json::to_string(&result["outcome"]["error"])?);Ok(())
}
