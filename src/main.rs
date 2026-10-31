use soil_fem_3d::{
    Error, Result, ground, mcc,
    mesh::Mesh,
    solver::{self, Case},
    stages,
};
use std::{env, fs, path::Path};

fn run() -> Result<()> {
    let args: Vec<_> = env::args().skip(1).collect();
    match args.as_slice() {
        [command, profile, output] if command == "spatial-convergence" => {
            if Path::new(output).exists() {
                return Err(Error::Input("spatial-convergence output already exists".into()));
            }
            let profile = match profile.as_str() {
                "ci" => soil_fem_3d::spatial_convergence::Profile::Ci,
                "user" => soil_fem_3d::spatial_convergence::Profile::User,
                _ => {
                    return Err(Error::Input(
                        "spatial-convergence profile must be ci or user".into(),
                    ));
                }
            };
            let report = soil_fem_3d::spatial_convergence::run(profile)?;
            write_new(output, &report)?;
            println!(
                "Spatial-convergence report saved. general_3d_demonstrated={} p6_complete={} smooth_elasticity={} uniform_srm_within_spacing={}",
                report.general_3d_spatial_convergence_demonstrated,
                report.p6_complete,
                report.smooth_linear_elasticity_orders_demonstrated,
                report.uniform_srm_factor_within_search_spacing
            );
        }
        [command,mesh,ground,case,spec,interval,output] if command=="prepare-mc-epoch-srm-gate-job" => {
            let interval:usize=interval.parse().map_err(|_|Error::Input("epoch SRM gate accepted-unit interval must be positive".into()))?;
            let m=soil_fem_3d::epoch_srm_gate_job::prepare(Path::new(mesh),Path::new(ground),Path::new(case),Path::new(spec),interval,Path::new(output))?;
            println!("Prepared unreduced epoch SRM stage-gate job {}",m.job_id);
        }
        [command,previous,output] if command=="prepare-mc-epoch-srm-gate-recovery" => {
            let m=soil_fem_3d::epoch_srm_gate_job::prepare_recovery(Path::new(previous),Path::new(output))?;
            println!("Prepared epoch SRM stage-gate recovery {}",m.job_id);
        }
        [command,folder] if command=="run-mc-epoch-srm-gate-job" => {
            match soil_fem_3d::epoch_srm_gate_job::run(Path::new(folder),None)? {
                soil_fem_3d::epoch_srm_gate_job::Outcome::Completed=>println!("Epoch SRM stage gate finished; read decision, gate effect and requirement. Job completion is not a safety certification"),
                soil_fem_3d::epoch_srm_gate_job::Outcome::Paused=>{println!("Epoch SRM stage gate paused at an accepted unit");std::process::exit(2);}
            }
        }
        [command,folder,budget] if command=="run-mc-epoch-srm-gate-job" => {
            let budget:usize=budget.parse().map_err(|_|Error::Input("epoch SRM gate advance budget must be positive".into()))?;
            match soil_fem_3d::epoch_srm_gate_job::run(Path::new(folder),Some(budget))? {
                soil_fem_3d::epoch_srm_gate_job::Outcome::Completed=>println!("Epoch SRM stage gate finished; read decision, gate effect and requirement. Job completion is not a safety certification"),
                soil_fem_3d::epoch_srm_gate_job::Outcome::Paused=>{println!("Epoch SRM stage gate paused at an accepted unit");std::process::exit(2);}
            }
        }
        [command,folder] if command=="stop-mc-epoch-srm-gate-job" => {
            println!("Stop requested: {}",soil_fem_3d::epoch_srm_gate_job::request_stop(Path::new(folder))?);
        }
        [command,folder] if command=="mc-epoch-srm-gate-job-status" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::epoch_srm_gate_job::status(Path::new(folder))?)?);
        }
        [command,folder] if command=="verify-mc-epoch-srm-gate-job" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::epoch_srm_gate_job::verify_snapshot(Path::new(folder))?)?);
        }
        [command,mesh,ground,case,spec,output] if command=="solve-mc-epoch-srm-gate" => {
            if Path::new(output).exists() {
                return Err(Error::Input("epoch SRM gate output already exists".into()));
            }
            let mesh: Mesh=serde_json::from_str(&fs::read_to_string(mesh)?)?;
            let ground: ground::Ground=serde_json::from_str(&fs::read_to_string(ground)?)?;
            let case: stages::StageCase=serde_json::from_str(&fs::read_to_string(case)?)?;
            let spec: soil_fem_3d::epoch_srm_gate::Specification=serde_json::from_str(&fs::read_to_string(spec)?)?;
            let report=soil_fem_3d::epoch_srm_gate::run(&mesh,&ground,&case,&spec)?;
            write_new(output,&report)?;
            println!("Epoch SRM stage gate saved; read decision, requirement and stop reason. This is not a safety certification");
        }
        [command,mesh,origin,spec,interval,output] if command=="prepare-mc-epoch-srm-job" => {
            let interval:usize=interval.parse().map_err(|_|Error::Input("epoch SRM material-attempt interval must be positive".into()))?;
            let m=soil_fem_3d::epoch_srm_job::prepare(Path::new(mesh),Path::new(origin),Path::new(spec),interval,Path::new(output))?;
            println!("Prepared accepted-origin epoch SRM job {}",m.job_id);
        }
        [command,previous,output] if command=="prepare-mc-epoch-srm-recovery" => {
            let m=soil_fem_3d::epoch_srm_job::prepare_recovery(Path::new(previous),Path::new(output))?;
            println!("Prepared epoch SRM recovery {}",m.job_id);
        }
        [command,folder] if command=="run-mc-epoch-srm-job" => {
            match soil_fem_3d::epoch_srm_job::run(Path::new(folder),None)? {
                soil_fem_3d::epoch_srm_job::Outcome::Completed=>println!("Epoch SRM finished; inspect physical evidence and requirement"),
                soil_fem_3d::epoch_srm_job::Outcome::Paused=>{println!("Epoch SRM paused at saved material-attempt boundary");std::process::exit(2);}
            }
        }
        [command,folder,budget] if command=="run-mc-epoch-srm-job" => {
            let budget:usize=budget.parse().map_err(|_|Error::Input("epoch SRM advance budget must be positive".into()))?;
            match soil_fem_3d::epoch_srm_job::run(Path::new(folder),Some(budget))? {
                soil_fem_3d::epoch_srm_job::Outcome::Completed=>println!("Epoch SRM finished; inspect physical evidence and requirement"),
                soil_fem_3d::epoch_srm_job::Outcome::Paused=>{println!("Epoch SRM paused at saved material-attempt boundary");std::process::exit(2);}
            }
        }
        [command,folder] if command=="stop-mc-epoch-srm-job" => {
            println!("Stop requested: {}",soil_fem_3d::epoch_srm_job::request_stop(Path::new(folder))?);
        }
        [command,folder] if command=="mc-epoch-srm-job-status" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::epoch_srm_job::status(Path::new(folder))?)?);
        }
        [command,folder] if command=="verify-mc-epoch-srm-job" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::epoch_srm_job::verify_snapshot(Path::new(folder))?)?);
        }
        [command,mesh,origin,spec,output] if command=="search-mc-epoch-srm" => {
            if Path::new(output).exists() {
                return Err(Error::Input("epoch SRM output already exists".into()));
            }
            let mesh: Mesh=serde_json::from_str(&fs::read_to_string(mesh)?)?;
            let origin: stages::Run=serde_json::from_str(&fs::read_to_string(origin)?)?;
            let spec: soil_fem_3d::epoch_srm_search::Specification=serde_json::from_str(&fs::read_to_string(spec)?)?;
            let report=soil_fem_3d::epoch_srm_search::run(&mesh,origin,spec)?;
            write_new(output,&report)?;
            println!("Accepted-origin epoch SRM saved; inspect requirement, stop reason and independent evidence");
        }
        [command,mesh,ground,case,spec,interval,output] if command=="prepare-srm-search-job" => {
            let interval:usize=interval.parse().map_err(|_|Error::Input("SRM accepted interval must be a positive integer".into()))?;
            let m=soil_fem_3d::srm_search_job::prepare(Path::new(mesh),Path::new(ground),Path::new(case),Path::new(spec),interval,Path::new(output))?;
            println!("Prepared immutable SRM search job {}",m.job_id);
        }
        [command,previous,output] if command=="prepare-srm-search-recovery" => {
            let m=soil_fem_3d::srm_search_job::prepare_recovery(Path::new(previous),Path::new(output))?;
            println!("Prepared SRM snapshot recovery {}",m.job_id);
        }
        [command,folder] if command=="run-srm-search-job" => {
            match soil_fem_3d::srm_search_job::run(Path::new(folder),None)? {
                soil_fem_3d::srm_search_job::Outcome::Completed=>println!("SRM search completed; inspect requirement and limitations separately"),
                soil_fem_3d::srm_search_job::Outcome::Paused=>{println!("SRM search paused with verified recovery snapshot");std::process::exit(2);}
            }
        }
        [command,folder,budget] if command=="run-srm-search-job" => {
            let budget:usize=budget.parse().map_err(|_|Error::Input("SRM advance budget must be positive".into()))?;
            match soil_fem_3d::srm_search_job::run(Path::new(folder),Some(budget))? {
                soil_fem_3d::srm_search_job::Outcome::Completed=>println!("SRM search completed; inspect requirement and limitations separately"),
                soil_fem_3d::srm_search_job::Outcome::Paused=>{println!("SRM search paused with verified recovery snapshot");std::process::exit(2);}
            }
        }
        [command,folder] if command=="stop-srm-search-job" => {
            println!("Cooperative SRM stop requested: {}",soil_fem_3d::srm_search_job::request_stop(Path::new(folder))?);
        }
        [command,folder] if command=="srm-search-job-status" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::srm_search_job::status(Path::new(folder))?)?);
        }
        [command,folder] if command=="verify-srm-search-job" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::srm_search_job::verify_snapshot(Path::new(folder))?)?);
        }
        [command,mesh,problem,output] if command=="mc-virtual-work-bound" => {
            if Path::new(output).exists() {
                return Err(Error::Input("virtual-work output already exists".into()));
            }
            let mesh: Mesh=serde_json::from_str(&fs::read_to_string(mesh)?)?;
            let problem: soil_fem_3d::srm_virtual_work::Problem=serde_json::from_str(&fs::read_to_string(problem)?)?;
            let report=soil_fem_3d::srm_virtual_work::audit(&mesh,&problem)?;
            write_new(output,&report)?;
            println!("Pointwise virtual-work assessment saved; no continuum stability certification inferred");
        }
        [command,mesh,ground,case,spec,output] if command=="search-mc-srm" => {
            if Path::new(output).exists() {
                return Err(Error::Input("SRM output already exists".into()));
            }
            let mesh: Mesh=serde_json::from_str(&fs::read_to_string(mesh)?)?;
            let ground: ground::Ground=serde_json::from_str(&fs::read_to_string(ground)?)?;
            let case: stages::StageCase=serde_json::from_str(&fs::read_to_string(case)?)?;
            let spec: soil_fem_3d::srm_search::Specification=serde_json::from_str(&fs::read_to_string(spec)?)?;
            let report=soil_fem_3d::srm_search::run(&mesh,ground,case,spec)?;
            write_new(output,&report)?;
            println!("Fresh-history SRM search saved; inspect physical interval, requirement and limitations");
        }
        [command,mesh,problem,output] if command=="mc-mean-stress-bound" => {
            let mesh: Mesh=serde_json::from_str(&fs::read_to_string(mesh)?)?;
            let problem: soil_fem_3d::srm_bound::Problem=serde_json::from_str(&fs::read_to_string(problem)?)?;
            let assessment=soil_fem_3d::srm_bound::analyze(&mesh,&problem)?;
            write_new(output,&assessment)?;
            println!("Mean-stress exclusion assessment saved; no collapse or stability classification inferred");
        }
        [command,mesh,origin,spec,output] if command=="prepare-mc-epoch-job" => {
            let m=soil_fem_3d::mc_epoch_job::prepare(Path::new(mesh),Path::new(origin),Path::new(spec),Path::new(output))?;
            println!("Prepared material continuation job {}",m.job_id);
        }
        [command,previous,output] if command=="prepare-mc-epoch-recovery" => {
            let m=soil_fem_3d::mc_epoch_job::prepare_recovery(Path::new(previous),Path::new(output))?;
            println!("Prepared material continuation recovery {}",m.job_id);
        }
        [command,folder] if command=="run-mc-epoch-job" => {
            mc_epoch_outcome(soil_fem_3d::mc_epoch_job::run(Path::new(folder),None)?)?;
        }
        [command,folder,budget] if command=="run-mc-epoch-job" => {
            let budget:usize=budget.parse().map_err(|_|Error::Input("MC attempt budget must be a positive integer".into()))?;
            mc_epoch_outcome(soil_fem_3d::mc_epoch_job::run(Path::new(folder),Some(budget))?)?;
        }
        [command,folder] if command=="stop-mc-epoch-job" => {
            println!("Cooperative material continuation stop requested: {}",soil_fem_3d::mc_epoch_job::request_stop(Path::new(folder))?);
        }
        [command,folder] if command=="mc-epoch-job-status" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::mc_epoch_job::status(Path::new(folder))?)?);
        }
        [command,folder] if command=="verify-mc-epoch-job" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::mc_epoch_job::verify_snapshot(Path::new(folder))?)?);
        }
        [command,mesh,ground,case,interval,output] if command=="prepare-durable-stages-job" => {
            let interval:usize=interval.parse().map_err(|_|Error::Input("checkpoint interval must be a positive integer".into()))?;
            let manifest=soil_fem_3d::stage_job::prepare_durable(Path::new(mesh),Path::new(ground),Path::new(case),interval,Path::new(output))?;
            println!("Prepared durable staged mechanics job {}",manifest.job_id);
        }
        [command,previous,output] if command=="prepare-stages-recovery" => {
            let manifest=soil_fem_3d::stage_job::prepare_recovery(Path::new(previous),Path::new(output))?;
            println!("Prepared accepted-snapshot recovery {}",manifest.job_id);
        }
        [command,mesh,ground,case,output] if command=="prepare-stages-job" => {
            let manifest=soil_fem_3d::stage_job::prepare(Path::new(mesh),Path::new(ground),Path::new(case),Path::new(output))?;
            println!("Prepared staged mechanics job {}",manifest.job_id);
        }
        [command,previous,output] if command=="prepare-stages-resume" => {
            let manifest=soil_fem_3d::stage_job::prepare_resume(Path::new(previous),Path::new(output))?;
            println!("Prepared staged mechanics restart {}",manifest.job_id);
        }
        [command,folder] if command=="run-stages-job" => {
            match soil_fem_3d::stage_job::run(Path::new(folder))? {
                soil_fem_3d::stage_job::Outcome::Succeeded=>println!("Staged mechanics job completed"),
                soil_fem_3d::stage_job::Outcome::Checkpointed=>{
                    println!("Staged mechanics job paused at an accepted increment; checkpoint saved");
                    std::io::Write::flush(&mut std::io::stdout())?;
                    std::process::exit(2);
                },
                soil_fem_3d::stage_job::Outcome::Failed=>return Err(Error::Solve("staged mechanics job failed; inspect status.json and incomplete.json".into())),
            }
        }
        [command,folder] if command=="stop-stages-job" => {
            println!("Cooperative stop requested: {}",soil_fem_3d::stage_job::request_stop(Path::new(folder))?);
        }
        [command,folder] if command=="stages-job-status" => {
            println!("{}",serde_json::to_string_pretty(&soil_fem_3d::stage_job::status(Path::new(folder))?)?);
        }
        [command,mesh_path,ground_path,case_path,count,out] if command=="checkpoint-one-way-construction" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way_construction::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let count:usize=count.parse().map_err(|_|Error::Input("checkpoint count must be positive".into()))?;
            let run=soil_fem_3d::one_way_construction::solve_until(&mesh,ground,case,count)?;
            let saved=soil_fem_3d::one_way_construction::Checkpoint::seal(run)?;write_new(out,&saved)?;
        }
        [command,mesh_path,ground_path,case_path,saved_path,out] if command=="resume-one-way-construction" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way_construction::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let saved:soil_fem_3d::one_way_construction::Checkpoint=serde_json::from_slice(&fs::read(saved_path)?)?;
            let run=saved.resume(&mesh,ground,case,None)?;write_new(out,&run)?;
            if !run.complete{return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"construction resumed incomplete".into())));}
        }
        [command,mesh_path,ground_path,case_path,out] if command=="solve-one-way-construction" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way_construction::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let run=soil_fem_3d::one_way_construction::solve(&mesh,ground,case)?;write_new(out,&run)?;
            if !run.complete{return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"construction incomplete".into())));}
            println!("construction complete: {} frames, {} removals",run.frames.len(),run.excavations.len());
        }
        [command,mesh_path,ground_path,case_path,out] if command=="one-way-construction-identity" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way_construction::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            write_new(out,&serde_json::json!({"mesh_sha256":soil_fem_3d::checkpoint::hash_json(&mesh)?,"ground_sha256":soil_fem_3d::checkpoint::hash_json(&ground)?,"case_sha256":soil_fem_3d::checkpoint::hash_json(&case)?}))?;
        }
        [command,mesh_path,ground_path,case_path,count,out] if command=="checkpoint-one-way" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let count:usize=count.parse().map_err(|_|Error::Input("checkpoint count must be positive".into()))?;
            let run=soil_fem_3d::one_way::solve_until(&mesh,ground,case,count)?;
            let saved=soil_fem_3d::one_way::Checkpoint::seal(run)?;write_new(out,&saved)?;
            println!("ONE_WAY checkpoint saved at {} s",saved.run.last_committed.hydraulic.time_s);
        }
        [command,mesh_path,ground_path,case_path,saved_path,out] if command=="resume-one-way" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let saved:soil_fem_3d::one_way::Checkpoint=serde_json::from_slice(&fs::read(saved_path)?)?;
            let run=saved.resume(&mesh,ground,case,None)?;write_new(out,&run)?;
            if !run.complete{return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"ONE_WAY resumed incomplete".into())));}
            println!("ONE_WAY resumed to {} s",run.last_committed.hydraulic.time_s);
        }
        [command,mesh_path,ground_path,case_path,out] if command=="solve-one-way" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let run=soil_fem_3d::one_way::solve(&mesh,ground,case)?;write_new(out,&run)?;
            if !run.complete{return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"ONE_WAY incomplete".into())));}
            println!("ONE_WAY complete: {} frames, {} s",run.frames.len(),run.last_committed.hydraulic.time_s);
        }
        [command,mesh_path,ground_path,case_path,out] if command=="one-way-identity" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:soil_fem_3d::one_way::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            write_new(out,&serde_json::json!({"mesh_sha256":soil_fem_3d::checkpoint::hash_json(&mesh)?,"ground_sha256":soil_fem_3d::checkpoint::hash_json(&ground)?,"case_sha256":soil_fem_3d::checkpoint::hash_json(&case)?}))?;
        }
        [command,mesh_path,case_path,out] if command=="richards-identity" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;mesh.validate()?;
            let case:soil_fem_3d::richards::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            write_new(out,&serde_json::json!({"mesh_sha256":soil_fem_3d::checkpoint::hash_json(&mesh)?,"case_sha256":soil_fem_3d::checkpoint::hash_json(&case)?}))?;
        }
        [command,mesh_path,case_path,count,out] if command=="checkpoint-richards" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::richards::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let count:usize=count.parse().map_err(|_|Error::Input("checkpoint count must be a positive integer".into()))?;
            let result=soil_fem_3d::richards::solve_until(&mesh,case,&vec![true;mesh.tetrahedra.len()],None,count)?;
            let checkpoint=soil_fem_3d::richards_checkpoint::Checkpoint::seal(result)?;
            write_new(out,&checkpoint)?;
            println!("Richards checkpoint saved at physical time {} s",checkpoint.run.last_committed.time_s);
        }
        [command,mesh_path,case_path,checkpoint_path,out] if command=="resume-richards" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::richards::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let checkpoint:soil_fem_3d::richards_checkpoint::Checkpoint=serde_json::from_slice(&fs::read(checkpoint_path)?)?;
            let result=checkpoint.resume(&mesh,case,&vec![true;mesh.tetrahedra.len()],None,None)?;
            write_new(out,&result)?;
            if !result.complete{return Err(Error::Solve(result.failure_reason.unwrap_or_else(||"resumed Richards flow incomplete".into())));}
            println!("Richards resumed to physical time {} s",result.last_committed.time_s);
        }
        [command,mesh_path,case_path,out] if command=="solve-richards" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::richards::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let result=soil_fem_3d::richards::solve(&mesh,case,&vec![true;mesh.tetrahedra.len()],None)?;
            write_new(out,&result)?;
            if !result.complete {return Err(Error::Solve(result.failure_reason.unwrap_or_else(||"Richards flow incomplete".into())));}
            println!("Richards flow completed: {} frames, {} attempts",result.frames.len(),result.attempts.len());
        }
        [command,mesh_path,case_path,out] if command=="solve-darcy" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::darcy::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let result=soil_fem_3d::darcy::solve(&mesh,&case,&vec![true;mesh.tetrahedra.len()],None)?;
            write_new(out,&result)?;
            if !result.complete {return Err(Error::Solve(result.failure_reason.unwrap_or_else(||"steady saturated flow incomplete".into())));}
            println!("Steady Darcy completed: {} material points, {} hydraulic regions",result.points.len(),result.balances.len());
        }
        [command,path,out] if command == "hydrostatic-face" => {
            let case:soil_fem_3d::hydrostatic::FaceCase=serde_json::from_slice(&fs::read(path)?)?;
            let result=soil_fem_3d::hydrostatic::run_face(case)?;
            write_new(out,&result)?;
            println!("Hydrostatic face integrated: {} cells, wet-domain error bound {:.3e} kN",result.visited_cells,result.domain_error_bound_kn);
        }
        [command,path,out] if command == "hardening-path" => {
            let path:soil_fem_3d::hardening::MaterialPath=serde_json::from_slice(&fs::read(path)?)?;
            let result=soil_fem_3d::hardening::run_path(path)?;
            write_new(out,&result)?;
            println!("IDH path completed: {} states",result.states.len());
        }
        [command,mesh_path,ground_path,case_path,spec_path,through,out] if command=="replay-through-stage"=>{
            let through:usize=through.parse().map_err(|_|Error::Input("stage gate count must be a nonnegative integer".into()))?;
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:stages::StageCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let spec:soil_fem_3d::replay::Specification=serde_json::from_slice(&fs::read(spec_path)?)?;
            let report=soil_fem_3d::replay::run_through_stage(&mesh,ground,case,spec,through)?;
            write_new(out,&report)?;
            if report.replay.trials.iter().any(|t|t.outcome!=soil_fem_3d::replay::Outcome::Completed) {return Err(Error::Solve("one or more stage-gate replays did not complete; no stability/Fs classification was inferred".into()));}
            println!("{} independent replays through{} construction stages completed; stability/Fs not evaluated",report.replay.trials.len(),through);
        }
        [command,mesh_path,ground_path,case_path,spec_path,out] if command=="replay-stages"=>{
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:stages::StageCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let spec:soil_fem_3d::replay::Specification=serde_json::from_slice(&fs::read(spec_path)?)?;
            let report=soil_fem_3d::replay::run(&mesh,ground,case,spec)?;
            write_new(out,&report)?;
            if report.trials.iter().any(|t|t.outcome!=soil_fem_3d::replay::Outcome::Completed) {return Err(Error::Solve("one or more replay trials did not complete; no stability/Fs classification was inferred".into()));}
            println!("{} independent construction replays completed; stability/Fs not evaluated",report.trials.len());
        }
        [command,mesh_path,out] if command=="mesh-identity"=>{
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;mesh.validate()?;
            write_new(out,&serde_json::json!({"schema":"soil-fem-mesh-identity/1","mesh_sha256":soil_fem_3d::checkpoint::hash_json(&mesh)?,"integration_points":mesh.integration_layout()?}))?;
        }
        [command,path,out] if command == "mc-path" => {
            let path:soil_fem_3d::mc::Path=serde_json::from_slice(&fs::read(path)?)?;
            let result=soil_fem_3d::mc::run_path(path)?;
            write_new(out,&result)?;
            println!("3D MC path completed: {} states",result.states.len());
        }
        [command,mesh_path,case_path,out] if command=="solve-mixed-mechanics" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::mechanics::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let result=soil_fem_3d::mixed_mechanics::solve(&mesh,case)?;
            write_new(out,&result)?;
            println!("Mixed mechanics completed: {} pressure DOFs, true free residual={:.3e}",result.pressure_dofs.len(),result.mechanics.true_relative_free_force_residual);
        }
        [command,mesh_path,case_path,out] if command=="solve-mechanics" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case:soil_fem_3d::mechanics::Case=serde_json::from_slice(&fs::read(case_path)?)?;
            let result=soil_fem_3d::mechanics::solve(&mesh,case)?;
            write_new(out,&result)?;
            println!("General mechanics completed: {} integration points, true free residual={:.3e}, constraint error={:.3e} m",result.integration_points.len(),result.true_relative_free_force_residual,result.maximum_constraint_error_m);
        }
        [command,mesh_path,ground_path,case_path,count,out] if command == "checkpoint-stages" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:stages::StageCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let count:usize=count.parse().map_err(|_|Error::Input("checkpoint count must be a positive integer".into()))?;
            let run=stages::solve_until(&mesh,ground,case,count)?;
            let checkpoint=soil_fem_3d::checkpoint::Checkpoint::seal(&mesh,run)?;
            write_new(out,&checkpoint)?;
            println!("Checkpoint saved at accepted stage {} progress {:.8}",checkpoint.run.last_committed.stage_index,checkpoint.run.last_committed.stage_progress);
        }
        [command,mesh_path,ground_path,case_path,checkpoint_path,out] if command == "resume-stages" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:stages::StageCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let checkpoint:soil_fem_3d::checkpoint::Checkpoint=serde_json::from_slice(&fs::read(checkpoint_path)?)?;
            let run=checkpoint.resume(&mesh,ground,case,None)?;
            write_new(out,&run)?;
            if !run.complete {return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"resumed run incomplete".into())));}
            println!("Resumed staged FEM completed: {} frames",run.frames.len());
        }
        [command,case_path,out] if command == "benchmark-tet4" => {
            let case:soil_fem_3d::benchmark::GridCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let report=soil_fem_3d::benchmark::tet4_grid(case)?;
            write_new(out,&report)?;
            println!("Tet4 linear benchmark: {} elements, {} free DOFs, {} CSR entries, {} GMRES iterations, true residual={:.3e}",report.tetrahedra,report.free_dofs,report.linear.stored_nonzeros,report.linear.iterations,report.linear.true_relative_residual);
        }
        [command,mesh_path,ground_path,case_path,out] if command == "solve-stages" => {
            let mesh:Mesh=serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground:ground::Ground=serde_json::from_slice(&fs::read(ground_path)?)?;
            let case:stages::StageCase=serde_json::from_slice(&fs::read(case_path)?)?;
            let run=stages::solve(&mesh,ground,case)?;
            write_new(out,&run)?;
            if !run.complete {return Err(Error::Solve(run.failure_reason.unwrap_or_else(||"staged run incomplete".into())));}
            println!("Soil stages completed: {} frames, {} accepted / {} rejected increments",run.frames.len(),run.attempts.iter().filter(|a|a.accepted).count(),run.attempts.iter().filter(|a|!a.accepted).count());
        }
        [command,mesh_path,ground_path,out] if command == "initialize-ground" => {
            let mesh: Mesh = serde_json::from_slice(&fs::read(mesh_path)?)?;
            let ground: ground::Ground = serde_json::from_slice(&fs::read(ground_path)?)?;
            let state = ground::initialize(&mesh,ground)?;
            write_new(out,&state)?;
            println!("initialized {} material points; interface triangles={}; true free residual={:.3e}",state.point_states.len(),state.shared_interface_triangles,state.true_relative_free_force_residual);
        }
        [command,path,out] if command == "mcc-path" => {
            let path: mcc::MaterialPath = serde_json::from_slice(&fs::read(path)?)?;
            let result = mcc::run_path(path)?;
            write_new(out,&result)?;
            println!("MCC path accepted: {} states",result.states.len());
        }
        [command,path] if command == "validate" => {
            let mesh: Mesh = serde_json::from_slice(&fs::read(path)?)?;
            let v = mesh.validate()?;
            println!("nodes={} tetrahedra={} volume_m3={v:.12}",mesh.nodes.len(),mesh.tetrahedra.len());
        }
        [command,mesh_path,case_path,out] if command == "solve-elastic" => {
            let mesh: Mesh = serde_json::from_slice(&fs::read(mesh_path)?)?;
            let case: Case = serde_json::from_slice(&fs::read(case_path)?)?;
            let result = solver::solve(&mesh,&case)?;
            write_new(out,&result)?;
            println!("iterations={} true_relative_residual={:.3e} equilibrium_error_kn={:?}",result.iterations,result.true_relative_residual,result.equilibrium_error_kn);
        }
        _ => return Err(Error::Input("usage: soil-fem-3d spatial-convergence ci|user NEW_REPORT.json | prepare-mc-epoch-srm-gate-job MESH.json GROUND.json CASE.json SPEC.json ACCEPTED_UNITS NEW_JOB | run-mc-epoch-srm-gate-job JOB [ADVANCE_BUDGET] | prepare-mc-epoch-srm-gate-recovery PREVIOUS_JOB NEW_JOB | stop-mc-epoch-srm-gate-job JOB | mc-epoch-srm-gate-job-status JOB | verify-mc-epoch-srm-gate-job JOB | solve-mc-epoch-srm-gate MESH.json GROUND.json CASE.json SPEC.json NEW_REPORT.json | prepare-mc-epoch-srm-job MESH.json ORIGINAL_RUN.json SPEC.json MATERIAL_ATTEMPTS_PER_SAVE NEW_JOB | run-mc-epoch-srm-job JOB [ADVANCE_BUDGET] | prepare-mc-epoch-srm-recovery PREVIOUS_JOB NEW_JOB | stop-mc-epoch-srm-job JOB | mc-epoch-srm-job-status JOB | verify-mc-epoch-srm-job JOB | search-mc-epoch-srm MESH.json ORIGINAL_RUN.json SPEC.json NEW_REPORT.json | prepare-srm-search-job MESH.json GROUND.json CASE.json SPEC.json ACCEPTED_INTERVAL NEW_JOB | run-srm-search-job JOB [ADVANCE_BUDGET] | prepare-srm-search-recovery PREVIOUS_JOB NEW_JOB | stop-srm-search-job JOB | srm-search-job-status JOB | verify-srm-search-job JOB | mc-virtual-work-bound MESH.json PROBLEM.json NEW_REPORT.json | search-mc-srm MESH.json GROUND.json CASE.json SPEC.json NEW_REPORT.json | mc-mean-stress-bound MESH.json PROBLEM.json NEW_REPORT.json | prepare-mc-epoch-job MESH.json ORIGINAL_RUN.json SPEC.json NEW_JOB | run-mc-epoch-job JOB [ATTEMPT_BUDGET] | prepare-mc-epoch-recovery PREVIOUS_JOB NEW_JOB | stop-mc-epoch-job JOB | mc-epoch-job-status JOB | verify-mc-epoch-job JOB | prepare-durable-stages-job MESH.json GROUND.json CASE.json ACCEPTED_INTERVAL NEW_JOB | prepare-stages-recovery PREVIOUS_JOB NEW_JOB | solve-mixed-mechanics MESH.json CASE.json NEW_RESULT.json | solve-one-way-construction MESH.json GROUND.json CASE.json NEW_RESULT.json | checkpoint-one-way-construction MESH.json GROUND.json CASE.json ACCEPTED_STEPS NEW_CHECKPOINT.json | resume-one-way-construction MESH.json GROUND.json CASE.json CHECKPOINT.json NEW_RESULT.json | solve-one-way MESH.json GROUND.json CASE.json NEW_RESULT.json | checkpoint-one-way MESH.json GROUND.json CASE.json ACCEPTED_STEPS NEW_CHECKPOINT.json | resume-one-way MESH.json GROUND.json CASE.json CHECKPOINT.json NEW_RESULT.json | solve-richards MESH.json CASE.json NEW_RESULT.json | checkpoint-richards MESH.json CASE.json ACCEPTED_STEPS NEW_CHECKPOINT.json | resume-richards MESH.json CASE.json CHECKPOINT.json NEW_RESULT.json | solve-darcy MESH.json FLOW.json NEW_RESULT.json | replay-through-stage MESH.json GROUND.json CASE.json SPEC.json THROUGH_STAGES NEW_REPORT.json | replay-stages MESH.json GROUND.json CASE.json SPEC.json NEW_REPORT.json | mesh-identity MESH.json NEW_IDENTITY.json | validate MESH.json | solve-elastic MESH.json CASE.json NEW_RESULT.json | solve-mechanics MESH.json CASE.json NEW_RESULT.json | initialize-ground MESH.json GROUND.json NEW_STATE.json | mc-path PATH.json NEW_RESULT.json | mcc-path PATH.json NEW_RESULT.json | hydrostatic-face CASE.json NEW_RESULT.json | hardening-path PATH.json NEW_RESULT.json | solve-stages MESH.json GROUND.json CASE.json NEW_RESULT.json | checkpoint-stages MESH.json GROUND.json CASE.json ACCEPTED_INCREMENTS NEW_CHECKPOINT.json | resume-stages MESH.json GROUND.json CASE.json CHECKPOINT.json NEW_RESULT.json | benchmark-tet4 CASE.json NEW_REPORT.json".into())),
    }
    Ok(())
}
fn mc_epoch_outcome(outcome: soil_fem_3d::mc_epoch_job::Outcome) -> Result<()> {
    use soil_fem_3d::mc_epoch_job::Outcome;
    match outcome {
        Outcome::Succeeded => println!("Material continuation completed; stability not evaluated"),
        Outcome::Paused => {
            println!("Material continuation paused; accepted epoch and next increment saved");
            std::io::Write::flush(&mut std::io::stdout())?;
            std::process::exit(2);
        }
        Outcome::Incomplete => return Err(Error::Solve(
            "material continuation incomplete; inspect report.json; no stability bound inferred"
                .into(),
        )),
        Outcome::Failed => {
            return Err(Error::Control(
                "material continuation failed; previous accepted recovery snapshot retained".into(),
            ));
        }
    }
    Ok(())
}
fn write_new(path: &str, value: &impl serde::Serialize) -> Result<()> {
    let path = Path::new(path);
    if path.exists() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            "output already exists",
        )
        .into());
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .ok_or_else(|| Error::Input("output must name a file".into()))?;
    let epoch = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|_| Error::Input("system clock precedes epoch".into()))?
        .as_nanos();
    let temporary = parent.join(format!(
        ".{}.{}.{}.tmp",
        name.to_string_lossy(),
        std::process::id(),
        epoch
    ));
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)?;
    let result = (|| -> Result<()> {
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, &mut file);
        serde_json::to_writer_pretty(&mut writer, value)?;
        std::io::Write::flush(&mut writer)?;
        drop(writer);
        file.sync_all()?;
        // Hard-link publication is atomic and fails if another process has
        // claimed the destination. Never rename over an existing result.
        fs::hard_link(&temporary, path)?;
        Ok(())
    })();
    drop(file);
    let _ = fs::remove_file(&temporary);
    result?;
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
