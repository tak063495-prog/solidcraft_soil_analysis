use soil_fem_3d::{Result,checkpoint::{self,Checkpoint},ground::{self,Ground},mesh::Mesh,
    geostatic::{self,PointSeed,Seed},stage_control::{Observer,Progress},stages::{self,StageCase}};
use serde_json::{Value,json};
use std::{cell::Cell,fs,path::Path};
struct Stop{queries:Cell<usize>,after:usize}
impl Stop{fn new(after:usize)->Self{Self{queries:Cell::new(0),after}}}
impl Observer for Stop{
    fn on_progress(&self,_:&Progress)->Result<()>{Ok(())}
    fn stop_requested(&self)->Result<bool>{let n=self.queries.get();self.queries.set(n+1);Ok(n>=self.after)}
}
fn scientific(mut v:Value)->Value{fn visit(v:&mut Value){match v{Value::Object(m)=>{m.remove("elapsed_seconds");for x in m.values_mut(){visit(x)}},Value::Array(a)=>for x in a{visit(x)},_=>()}}visit(&mut v);v}
fn save(path:&Path,v:&Value)->Result<()>{use std::io::Write;let mut f=fs::OpenOptions::new().write(true).create_new(true).open(path)?;f.write_all(&serde_json::to_vec_pretty(v)?)?;f.sync_all()?;Ok(())}
fn main()->Result<()>{
    let args:Vec<_>=std::env::args().collect();let runtime=Path::new(&args[1]);let out=Path::new(&args[2]);let fixed=args[3]=="fixed";
    fs::create_dir(out)?;let mut rows=vec![];
    for projected in [false,true]{
        let mut mesh:Mesh=serde_json::from_slice(&fs::read(runtime.join("examples/hardening-column-tet10-mesh.json"))?)?;
        let cut=mesh.tetrahedra.iter().position(|t|t.iter().filter(|&&i|mesh.nodes[i][2]==1.).count()==3).unwrap();
        mesh.excavation_sets.insert("LocalCut".into(),vec![cut]);
        let mut g:Ground=serde_json::from_slice(&fs::read(runtime.join("examples/ground-mixed-mcc.json"))?)?;
        for material in &mut g.materials{material.unit_weight_kn_m3=1e-8;
            let soil_fem_3d::soil::Model::Mcc(par)=&mut material.model else{panic!("MCC required")};
            par.critical_state_m=1.2;par.lambda=0.2;par.kappa=0.04;par.poisson=0.3;
            par.reference_pressure_kpa=1.;par.specific_volume_ncl_at_reference=3.5;}
        for layer in &mut g.layers{layer.k0_x=1.;layer.k0_y=1.;layer.consolidation=ground::Consolidation::MccPcOverMeanEffective{ratio:1.};}
        g.excavation_bodies=vec!["LocalCut".into()];
        let initial=ground::initialize(&mesh,g.clone())?;
        let points=initial.point_states.into_iter().map(|q|PointSeed{element_index:q.element_index,local_index:q.local_index,
            material_id:q.material_id,position_m:q.centroid_m,state:q.constitutive_state.into_state()}).collect();
        g.initialization.method="equilibrate".into();
        g.initialization.geostatic=Some(geostatic::Configuration{seed:Seed::IntegrationPoints{mesh_sha256:checkpoint::hash_json(&mesh)?,points},
            initial_boundary_factor:0.,maximum_initial_engineering_strain:0.05,relative_volume_tolerance:1e-8,
            base_loads:Default::default(),hydrostatic_face_sets:vec![],verification_reject_once:true});
        g.layers.clear();g.bodies=mesh.regions.iter().map(|r|geostatic::Body{body:r.body.clone(),material_id:r.material_id}).collect();
        let mut raw:Value=serde_json::from_slice(&fs::read(runtime.join("examples/stages-mixed-mcc.json"))?)?;
        raw["schema"]=json!(if projected{"soil-fem-projected-mcc-stages/1"}else{"soil-fem-mixed-mcc-stages/1"});
        raw["settings"]["initial_increment"]=json!(0.25);raw["settings"]["maximum_free_dofs"]=json!(2000);
        raw["settings"]["relative_force_tolerance"]=json!(1e-9);
        raw["settings"]["relative_linear_tolerance"]=json!(1e-11);
        raw["boundary"]["constraints"][3]["value_m"]=json!(-0.005);
        raw["stages"]=json!([{"kind":"boundary","name":"load","target_factor":1.},
            {"kind":"boundary","name":"unload","target_factor":0.8},
            {"kind":"excavate","name":"cut","bodies":["LocalCut"],"remaining_surface_pressure_kpa":0.},
            {"kind":"boundary","name":"after_cut_unload","target_factor":0.6}]);
        raw["pressure_on_original_surface_only"]=json!(true);raw["verification_reject_once_stages"]=json!(["cut"]);
        let case:StageCase=serde_json::from_value(raw)?;let dir=out.join(if projected{"projected"}else{"mixed"});fs::create_dir(&dir)?;
        for (name,v) in [("mesh.json",serde_json::to_value(&mesh)?),("ground.json",serde_json::to_value(&g)?),("case.json",serde_json::to_value(&case)?)]{save(&dir.join(name),&v)?;}
        let full=stages::solve(&mesh,g.clone(),case.clone())?;assert!(full.complete,"{:?}",full.failure_reason);
        save(&dir.join("uninterrupted.json"),&serde_json::to_value(&full)?)?;
        assert!(full.attempts.iter().any(|a|a.stage=="cut"&&!a.accepted&&a.committed_state_unchanged_on_rejection));
        let geo=full.geostatic.as_ref().unwrap();assert!(geo.complete);
        let n=geo.attempts.iter().filter(|a|a.accepted).count()+full.attempts.iter().filter(|a|a.accepted).position(|a|a.stage=="cut").unwrap()+1;
        let partial=stages::solve_until(&mesh,g.clone(),case.clone(),n)?;
        assert!(partial.last_committed.stage_progress>0.&&partial.last_committed.stage_progress<1.);
        assert_eq!(partial.last_committed.active_elements.iter().filter(|&&a|!a).count(),1);
        assert!(partial.last_committed.material_states.iter().any(|s|s.plastic_strain().iter().any(|v|v.abs()>1e-7)));
        assert!(!partial.excavation_audits.is_empty());
        assert!(partial.last_committed.excavation_release_load_kn.iter().any(|v|v.abs()>1e-8));
        let before=serde_json::to_vec(&partial)?;save(&dir.join("before.json"),&serde_json::to_value(&partial)?)?;
        let clone_run=||->Result<stages::Run>{Ok(serde_json::from_slice(&before)?)};
        let stopped=Checkpoint::seal(&mesh,clone_run()?)?.resume_observed(&mesh,g.clone(),case.clone(),&Stop::new(0))?;
        save(&dir.join("stopped.json"),&serde_json::to_value(&stopped)?)?;
        assert_eq!(serde_json::to_vec(&stopped)?==before,fixed);
        if fixed{
            let resumed=Checkpoint::seal(&mesh,stopped)?.resume(&mesh,g.clone(),case.clone(),None)?;
            assert!(resumed.complete);assert_eq!(scientific(serde_json::to_value(&resumed)?),scientific(serde_json::to_value(&full)?));
            save(&dir.join("resumed.json"),&serde_json::to_value(&resumed)?)?;
            let late=Checkpoint::seal(&mesh,clone_run()?)?.resume_observed(&mesh,g.clone(),case.clone(),&Stop::new(1))?;
            assert_eq!(late.failure_reason.as_deref(),Some("checkpoint_requested"));
            assert_eq!(late.last_committed.stage_index,partial.last_committed.stage_index);
            assert!(late.last_committed.stage_progress>partial.last_committed.stage_progress);
            assert_eq!(late.last_committed.active_elements,partial.last_committed.active_elements);
            for point in mesh.integration_layout()?.into_iter().enumerate().filter_map(|(i,q)|(q.element_index==cut).then_some(i)){
                assert_eq!(serde_json::to_vec(&late.last_committed.material_states[point])?,serde_json::to_vec(&partial.last_committed.material_states[point])?);
            }
            save(&dir.join("late-stop.json"),&serde_json::to_value(&late)?)?;
            let resumed=Checkpoint::seal(&mesh,late)?.resume(&mesh,g.clone(),case.clone(),None)?;
            assert!(resumed.complete);assert_eq!(scientific(serde_json::to_value(&resumed)?),scientific(serde_json::to_value(&full)?));
            save(&dir.join("late-stop-resumed.json"),&serde_json::to_value(&resumed)?)?;
        }else{
            assert_eq!(stopped.last_committed.stage_index,0);
            assert!(stopped.last_committed.active_elements.iter().all(|&a|a));
            assert!(stopped.excavation_audits.is_empty());assert!(Checkpoint::seal(&mesh,stopped).is_err());
        }
        rows.push(json!({"projected":projected,"fixed":fixed,"actual_local_removed_elements":1,
            "plastic_and_release_and_rejected_cut_attempt_present":true,"saved_stage_progress":partial.last_committed.stage_progress,
            "stopped_all_Run_bytes_equal":fixed,"early_and_late_resumed_entire_history_equal_except_elapsed":fixed}));
    }
    save(&out.join("terminal.json"),&json!({"passed":true,"fixed":fixed,"rows":rows,"p6_complete":false}))?;Ok(())
}
