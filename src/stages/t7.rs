//! Isolated evaluations for continuation on the dry A-01 boundary stage.
//! The caller owns the accepted checkpoint. An evaluation cannot commit it.
use super::*;
use serde_json::{Value, json};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Request {
    parent: Committed,
    stage: usize,
    origin_factor: f64,
    progress: f64,
    free_m: Option<Vec<f64>>,
    contact_subdivision: usize,
}

pub fn evaluate(mesh: &Mesh, ground: &Ground, case: &StageCase, request: Request) -> Result<Value> {
    let model = Model::new(mesh, ground, case)?;
    let Request {
        parent,
        stage,
        origin_factor,
        progress,
        free_m,
        contact_subdivision,
    } = request;
    let Some(Stage::Boundary { target_factor, .. }) =
        stage.checked_sub(1).and_then(|s| case.stages.get(s))
    else {
        return Err(Error::Input(
            "T7 requires a proportional boundary stage".into(),
        ));
    };
    if !origin_factor.is_finite()
        || !progress.is_finite()
        || !(1..=64).contains(&contact_subdivision)
        || mesh.tetrahedra.len() > 144
        || model.beams.is_some()
        || model.anchors.is_some()
        || model.wet_joints.is_some()
        || parent.joints.is_none()
        || parent.pore_pressure_kpa.is_some()
        || !parent.active_elements.iter().all(|v| *v)
        || !(parent.stage_index == stage
            || (parent.stage_index.checked_add(1) == Some(stage) && parent.stage_progress == 1.))
    {
        return Err(Error::Input(
            "T7 dry parallel-interface continuation scope".into(),
        ));
    }
    validate_committed(&model, &parent, case)?;
    let settings = &case.settings;
    let factor = origin_factor + progress * (target_factor - origin_factor);
    let active = &parent.active_elements;
    let dofs = model.dofs(active, settings, factor, parent.joints.as_ref())?;
    let transform = dofs
        .transform
        .as_ref()
        .ok_or_else(|| Error::Input("T7 requires explicit supports".into()))?;
    let flat = parent
        .displacement_m
        .iter()
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    let free = free_m.unwrap_or_else(|| transform.reduce(&flat));
    if free.len() != dofs.count || !free.iter().all(|v| v.is_finite()) {
        return Err(Error::Input("T7 free coordinates".into()));
    }
    let flat: Vec<_> = transform
        .expand(&free)
        .iter()
        .zip(&transform.offset)
        .map(|(a, b)| a + b)
        .collect();
    let u: Vec<[f64; 3]> = flat.as_chunks::<3>().0.to_vec();
    let mut counts = model.predictor_partitions(&parent, &u, active, settings)?;
    for n in &mut counts[model.point_count..] {
        *n = n.saturating_mul(contact_subdivision);
    }
    if counts.iter().any(|&n| n > 50000) {
        return Err(Error::Solve("T7 integration partition budget".into()));
    }
    let load = model.physical_load(
        stage,
        active,
        parent.surface_pressure_kpa,
        case.pressure_on_original_surface_only,
        factor,
        &parent.element_material_ids,
        None,
        None,
        None,
    )?;
    let mut eval = model.evaluate(
        &parent,
        &u,
        active,
        &dofs,
        settings,
        EvaluationMode {
            tangent: true,
            factor: Some(factor),
            partitions: Some(&counts),
            ..EvaluationMode::plain()
        },
    )?;
    let geometry = crate::rigidity::Geometry::new(mesh, active)?;
    let support = model
        .joint_support_space(&geometry, active, &dofs, &eval)?
        .full();
    let reaction: Vec<_> = eval
        .internal
        .iter()
        .zip(&load)
        .map(|(a, b)| a - b)
        .collect();
    let residual = model.reduce(&reaction, &dofs);
    let tangent = eval
        .tangent
        .take()
        .ok_or_else(|| Error::Solve("missing trial tangent".into()))?;
    let columns = (0..dofs.count)
        .map(|j| {
            let mut unit = DVector::zeros(dofs.count);
            unit[j] = 1.;
            tangent.action(&unit).map(|c| c.as_slice().to_vec())
        })
        .collect::<Result<Vec<_>>>()?;
    let constraint_error = transform.maximum_constraint_error(&flat);
    let candidate = Committed {
        displacement_m: u,
        material_states: eval.states,
        joints: eval.joints.map(|j| j.state),
        applied_nodal_load_kn: load,
        boundary_factor: factor,
        stage_index: stage,
        stage_progress: progress,
        ..parent
    };
    Ok(
        json!({"candidate":candidate,"residual_kn":residual.as_slice(),"free_m":free,"dof_rows":dofs.rows,
              "offset_m":transform.offset,"tangent_columns":columns,"partitions":counts,"reaction_kn":reaction,
              "supported":support,"constraint_error_m":constraint_error}),
    )
}
