//! Read-only force and tangent actions at the committed origin of an increment.
use super::*;

#[derive(Serialize)]
pub struct CheckpointResponse {
    dof_rows: Vec<Vec<(usize, f64)>>,
    reduced_residual_kn: Vec<f64>,
    actions_kn: Vec<Vec<f64>>,
}

impl Run {
    /// Evaluate the production partitioned material map at zero strain increment.
    /// Columns are supplied in reduced displacement coordinates. No trial is solved.
    pub fn checkpoint_response(
        &self,
        mesh: &Mesh,
        directions: &[Vec<f64>],
    ) -> Result<CheckpointResponse> {
        mesh.validate()?;
        self.ground.validate()?;
        self.case.validate()?;
        if crate::checkpoint::hash_json(mesh)? != self.provenance.mesh_sha256 {
            return Err(Error::Input(
                "mesh differs from the recorded checkpoint".into(),
            ));
        }
        let mut model = Model::new(mesh, &self.ground, &self.case)?;
        model
            .reference_displacement_m
            .clone_from(&self.reference_displacement_m);
        let state = &self.last_committed;
        validate_committed(&model, state, &self.case)?;
        let settings = &self.case.settings;
        let dofs = model.dofs(
            &state.active_elements,
            settings,
            state.boundary_factor,
            state.joints.as_ref(),
        )?;
        if directions.is_empty()
            || directions
                .iter()
                .any(|v| v.len() != dofs.count || !v.iter().all(|x| x.is_finite()))
        {
            return Err(Error::Input(
                "expected finite reduced displacement directions".into(),
            ));
        }
        let partitions =
            vec![1; model.point_count + model.joints.as_ref().map_or(0, |j| j.point_count())];
        let eval = model.evaluate(
            state,
            &state.displacement_m,
            &state.active_elements,
            &dofs,
            settings,
            EvaluationMode {
                tangent: true,
                partitions: Some(&partitions),
                factor: Some(state.boundary_factor),
                ..EvaluationMode::plain()
            },
        )?;
        let residual = state
            .applied_nodal_load_kn
            .iter()
            .zip(eval.internal)
            .map(|(f, i)| f - i)
            .collect::<Vec<_>>();
        let reduced = model.reduce(&residual, &dofs);
        let tangent = eval
            .tangent
            .ok_or_else(|| Error::Solve("missing checkpoint tangent".into()))?;
        let actions_kn = directions
            .iter()
            .map(|v| {
                tangent
                    .action(&DVector::from_column_slice(v))
                    .map(|y| y.as_slice().to_vec())
            })
            .collect::<Result<_>>()?;
        Ok(CheckpointResponse {
            dof_rows: dofs.rows,
            reduced_residual_kn: reduced.as_slice().to_vec(),
            actions_kn,
        })
    }
}
