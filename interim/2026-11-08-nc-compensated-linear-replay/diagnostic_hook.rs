// Appended only to an isolated copy of stages.rs. No continuation is produced.
// A single recorded increment is evaluated; this is not checkpoint resume.
pub fn diagnose_recorded_nc_trial(
    mesh: &Mesh,
    saved: &crate::checkpoint::Checkpoint,
    next: f64,
    observer: &dyn crate::stage_control::Observer,
) -> Result<serde_json::Value> {
    let run = &saved.run;
    let case = &run.case;
    let ground = &run.ground;
    if saved.schema != "soil-fem-checkpoint/1"
        || saved.run_sha256 != crate::checkpoint::hash_json(run)?
        || saved.identity != run.provenance
        || saved.identity.mesh_sha256 != crate::checkpoint::hash_json(mesh)?
        || saved.identity.ground_sha256 != crate::checkpoint::hash_json(ground)?
        || saved.identity.case_sha256 != crate::checkpoint::hash_json(case)?
        || saved.identity.solver_binary_sha256
            != "598a0d42a13bb00bdaffef52e82c45debb2a974c0f54c94cce9926b27fc4bfd6"
    {
        return Err(Error::Input("diagnostic input/history identity mismatch".into()));
    }
    let old = &run.last_committed;
    let c = run.continuation.as_ref().ok_or_else(|| Error::Input("missing continuation".into()))?;
    if mesh.nodes.len() != 54 || mesh.tetrahedra.len() != 12
        || old.stage_index != 3 || c.stage_index != 2
        || !(0.9513997..=0.9514001).contains(&old.stage_progress)
        || !(next > old.stage_progress && next <= 0.951401)
        || case.settings.relative_force_tolerance != 1e-8
        || case.settings.relative_linear_tolerance != 1e-10
        || case.settings.mcc_integration != Some(crate::mcc::Integration::Midpoint)
        || run.complete
    {
        return Err(Error::Input("diagnostic restricted to recorded NC stop window".into()));
    }
    let before = serde_json::to_vec(saved)?;
    let mut model = Model::new(mesh, ground, case)?;
    model.observer = Some(observer);
    model.reference_displacement_m = run.reference_displacement_m.clone();
    let target = Target {
        joint_water: old.joint_water.clone(),
        pore_pressure_kpa: old.pore_pressure_kpa.clone(),
        water_level: old.hydrostatic_water_level_m,
        anchor_control: None,
        active: c.active.clone(),
        load: c.end_load.iter().zip(&c.holding_load).map(|(f,h)| f+(1.-next)*h).collect(),
        release: c.holding_load.iter().map(|h| (1.-next)*h).collect(),
        pressure: old.surface_pressure_kpa,
        boundary_factor: old.boundary_factor,
        stage: 3,
        progress: next,
    };
    let outcome = match model.trial(old, target, &case.settings) {
        Ok(t) => serde_json::json!({"converged": true, "state": t.state,
            "iterations": t.iterations, "relative_force_residual": t.residual}),
        Err(e) => serde_json::json!({"converged": false, "error": e.to_string()}),
    };
    if before != serde_json::to_vec(saved)? {
        return Err(Error::Control("diagnostic mutated saved history".into()));
    }
    Ok(serde_json::json!({"passed": false, "diagnostic_only": true,
        "checkpoint_resume": false, "new_continuation_created": false,
        "source_progress": old.stage_progress, "target_progress": next,
        "saved_checkpoint_unchanged": true, "outcome": outcome}))
}
