//! Read-only reconstruction of an accepted dry construction endpoint for SRM.
use super::*;

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Endpoint {
    pub physical_force_kn: Vec<f64>,
    pub physical_force_sha256: String,
    pub active_elements: Vec<bool>,
    pub boundary_constraints: Vec<crate::boundary::Constraint>,
    pub maximum_relative_reconstructed_load_difference: f64,
    pub audit: Frame,
}

pub fn inspect(mesh: &Mesh, run: &Run) -> Result<Endpoint> {
    // Reuse the exact existing origin and full committed-state validation.
    let _validated = McParameterMechanics::new(mesh, run)?;
    let mut model = Model::new(mesh, &run.ground, &run.case)?;
    model.reference_displacement_m = run.reference_displacement_m.clone();
    let state = &run.last_committed;
    if state
        .pore_pressure_kpa
        .as_ref()
        .is_some_and(|p| p.iter().any(|v| *v != 0.))
        || state.hydrostatic_water_level_m.is_some()
        || state
            .element_material_ids
            .iter()
            .enumerate()
            .any(|(i, id)| *id != mesh.material_ids[i])
    {
        return Err(Error::Input(
            "SRM endpoint requires unchanged MC material IDs and zero water pressure".into(),
        ));
    }
    let physical = model.physical_load(
        &state.active_elements,
        state.surface_pressure_kpa,
        run.case.pressure_on_original_surface_only,
        state.boundary_factor,
        &state.element_material_ids,
        None,
        None,
        None,
    )?;
    let scale = physical.iter().map(|x| x.abs()).fold(1., f64::max);
    let difference = physical
        .iter()
        .zip(&state.applied_nodal_load_kn)
        .map(|(a, b)| (a - b).abs() / scale)
        .fold(0., f64::max);
    if !difference.is_finite() || difference > 1e-12 {
        return Err(Error::Input("accepted endpoint load differs from reconstructed physical load; residual holding/release loads cannot define an SRM witness".into()));
    }
    let audit = model.frame(
        "independent_srm_endpoint_audit".into(),
        state,
        state,
        &run.case,
    )?;
    if !audit.relative_force_residual.is_finite()
        || audit.relative_force_residual > run.case.settings.relative_force_tolerance
        || !audit.maximum_normalized_yield.is_finite()
        || audit.maximum_normalized_yield > 1e-8
        || !audit.maximum_constraint_error_m.is_finite()
        || audit.maximum_constraint_error_m > 1e-10
    {
        return Err(Error::Input(
            "SRM endpoint failed reconstructed force/yield/constraint audit".into(),
        ));
    }
    let constraints = run
        .case
        .boundary
        .as_ref()
        .ok_or_else(|| Error::Input("SRM requires explicit boundary constraints".into()))?
        .constraints
        .clone();
    if constraints
        .iter()
        .any(|c| c.value_m != 0. || c.gradient != [0.; 3])
    {
        return Err(Error::Input(
            "this SRM witness requires zero prescribed displacements".into(),
        ));
    }
    Ok(Endpoint {
        physical_force_sha256: crate::checkpoint::hash_json(&physical)?,
        physical_force_kn: physical,
        active_elements: state.active_elements.clone(),
        boundary_constraints: constraints,
        maximum_relative_reconstructed_load_difference: difference,
        audit,
    })
}

/// Re-audit an accepted parameter epoch against the unchanged construction
/// origin. This never converts an epoch into a normal staged Run/Checkpoint.
pub fn inspect_epoch(mesh: &Mesh, origin: &Run, saved: &McEpochCheckpoint) -> Result<Endpoint> {
    let validated = McParameterMechanics::resume(mesh, origin, saved)?;
    let mut endpoint = inspect(mesh, origin)?;
    let mut model = Model::new(mesh, &origin.ground, &origin.case)?;
    model.reference_displacement_m = origin.reference_displacement_m.clone();
    model.material_models = validated
        .parameters()
        .iter()
        .cloned()
        .map(crate::soil::Model::MohrCoulomb)
        .collect();
    // resume already validates these frozen metadata and all original force,
    // yield and constraint limits. Keep a complete current-material Frame.
    endpoint.audit = model.frame(
        "independent_srm_epoch_endpoint_audit".into(),
        validated.committed(),
        &origin.last_committed,
        &origin.case,
    )?;
    Ok(endpoint)
}
