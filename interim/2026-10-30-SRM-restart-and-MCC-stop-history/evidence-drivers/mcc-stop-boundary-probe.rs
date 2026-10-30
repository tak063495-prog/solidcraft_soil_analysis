use soil_fem_3d::{Result, checkpoint::Checkpoint, ground::Ground, mesh::Mesh,
    stage_control::{Observer,Progress}, stages::{self,StageCase}};
use serde_json::{Value,json};
use std::{cell::Cell,fs,path::Path};
struct Stop { queries:Cell<usize>, events:Cell<usize>, after:usize }
impl Stop { fn new(after:usize)->Self {Self{queries:Cell::new(0),events:Cell::new(0),after}} }
impl Observer for Stop {
    fn on_progress(&self,_:&Progress)->Result<()> {self.events.set(self.events.get()+1);Ok(())}
    fn stop_requested(&self)->Result<bool>{let n=self.queries.get();self.queries.set(n+1);Ok(n>=self.after)}
}
fn scientific(mut v:Value)->Value {
    fn visit(v:&mut Value){match v {Value::Object(m)=>{m.remove("elapsed_seconds");for x in m.values_mut(){visit(x)}},
        Value::Array(a)=>for x in a{visit(x)},_=>()}}
    visit(&mut v);v
}
fn save(path:&Path,v:&Value)->Result<()> {use std::io::Write;let mut f=fs::OpenOptions::new().write(true).create_new(true).open(path)?;
    f.write_all(&serde_json::to_vec_pretty(v)?)?;f.sync_all()?;Ok(())}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args().collect();let source=Path::new(&args[1]);let out=Path::new(&args[2]);
    let fixed=args[3]=="fixed";fs::create_dir(out)?;let mut rows=vec![];
    for mixed in [false,true] {
        let mut mesh:Mesh=serde_json::from_slice(&fs::read(source.join("examples/hardening-column-tet10-mesh.json"))?)?;
        mesh.face_sets.insert("top".into(),soil_fem_3d::boundary::surface_owners(&mesh).into_iter()
            .filter(|(face,owners)|owners.len()==1&&face.iter().all(|&i|mesh.nodes[i][2]==1.))
            .map(|(face,_)|face).collect());
        let g:Ground=serde_json::from_slice(&fs::read(source.join("examples/ground-geostatic-mcc-column.json"))?)?;
        let mut c:StageCase=serde_json::from_slice(&fs::read(source.join("examples/stages-geostatic-column.json"))?)?;
        if mixed {c.schema="soil-fem-mixed-mcc-stages/1".into()}
        c.settings.initial_increment=0.25;
        c.verification_reject_once_stages=vec!["small_reload".into()];
        let dir=out.join(if mixed {"mixed"} else {"displacement"});fs::create_dir(&dir)?;
        save(&dir.join("mesh.json"),&serde_json::to_value(&mesh)?)?;
        save(&dir.join("ground.json"),&serde_json::to_value(&g)?)?;
        save(&dir.join("case.json"),&serde_json::to_value(&c)?)?;
        let full=stages::solve(&mesh,g.clone(),c.clone())?;assert!(full.complete,"{:?}",full.failure_reason);
        let geo=full.geostatic.as_ref().unwrap();assert!(geo.complete);
        assert!(geo.last_committed.material_states.iter().any(|p|p.as_mcc().unwrap().preconsolidation_pressure_kpa>140.01));
        assert!(full.attempts.iter().any(|a|!a.accepted&&a.rejected_trial_plastic_changed&&a.committed_state_unchanged_on_rejection));
        let count=geo.attempts.iter().filter(|a|a.accepted).count();
        let mut cases=vec![];
        for (name,n) in [("completed-geostatic",count),("plastic-construction",count+1)] {
            let saved=stages::solve_until(&mesh,g.clone(),c.clone(),n)?;
            assert_eq!(saved.failure_reason.as_deref(),Some("checkpoint_requested"));
            if n>count {assert!(saved.last_committed.stage_index>0&&saved.last_committed.stage_progress<1.);
                assert!(saved.last_committed.material_states.iter().any(|p|p.plastic_strain().iter().any(|x|x.abs()>1e-8)));}
            let before=serde_json::to_value(&saved)?;save(&dir.join(format!("{name}-before.json")),&before)?;
            let trace=Stop::new(0);
            let stopped=Checkpoint::seal(&mesh,saved)?.resume_observed(&mesh,g.clone(),c.clone(),&trace)?;
            let exact=serde_json::to_value(&stopped)?==before;
            save(&dir.join(format!("{name}-stopped.json")),&serde_json::to_value(&stopped)?)?;
            assert_eq!(trace.events.get(),0);
            assert_eq!(exact,fixed,"the old boundary fault must be reproduced, and fixed candidate must keep every byte");
            if fixed {
                let checkpoint=Checkpoint::seal(&mesh,stopped)?;
                let bytes=serde_json::to_vec(&checkpoint)?;
                let checkpoint:Checkpoint=serde_json::from_slice(&bytes)?;
                let resumed=checkpoint.resume(&mesh,g.clone(),c.clone(),None)?;
                assert!(resumed.complete);
                assert_eq!(scientific(serde_json::to_value(&resumed)?),scientific(serde_json::to_value(&full)?));
                save(&dir.join(format!("{name}-resumed.json")),&serde_json::to_value(&resumed)?)?;
            } else {
                assert_eq!(stopped.failure_reason.as_deref(),Some("geostatic initialization incomplete"));
                assert_eq!(stopped.last_committed.stage_index,0);
                assert!(Checkpoint::seal(&mesh,stopped).is_err());
            }
            cases.push(json!({"boundary":name,"stopped_run_byte_exact":exact,
                "all_scientific_history_matches_uninterrupted_if_fixed":fixed,"events_after_stop":trace.events.get()}));
        }
        save(&dir.join("uninterrupted.json"),&serde_json::to_value(&full)?)?;
        rows.push(json!({"mixed":mixed,"elements":mesh.tetrahedra.len(),"integration_points":full.last_committed.material_states.len(),
            "actual_geostatic_hardening_and_rejected_plastic_trial":true,"cases":cases}));
    }
    save(&out.join("terminal.json"),&json!({"passed":true,"fixed_candidate":fixed,"rows":rows,
        "only_elapsed_seconds_removed_from_final_scientific_comparison":true,"p6_complete":false}))?;
    Ok(())
}
