//! Opt-in step doubling of GLOBAL equilibria. Probe states are never committed.
use super::*;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathErrorAudit {
    pub midpoint_progress: f64,
    pub base_increment: f64,
    pub actual_increment: f64,
    pub slip_budget_m: f64,
    pub dissipation_budget_kpa_m: f64,
    pub maximum_slip_error_m: f64,
    pub maximum_dissipation_error_kpa_m: f64,
    pub slip_errors_m: Vec<[f64; 3]>,
    pub dissipation_errors_kpa_m: Vec<f64>,
    pub old_status: Vec<crate::joint_law::Status>,
    pub midpoint_status: Vec<crate::joint_law::Status>,
    pub coarse_status: Vec<crate::joint_law::Status>,
    pub fine_status: Vec<crate::joint_law::Status>,
    pub hidden_midpoint_switch_points: usize,
    pub endpoint_disagreement_points: usize,
    pub midpoint_relative_force_residual: f64,
    pub fine_relative_force_residual: f64,
    pub midpoint_linear_solves: Vec<LinearReport>,
    pub fine_linear_solves: Vec<LinearReport>,
    pub reject: bool,
}

impl Model<'_> {
    pub(super) fn trial_with_path_error(
        &self,
        old: &Committed,
        target: Target,
        settings: &Settings,
    ) -> Result<(Trial, Option<PathErrorAudit>)> {
        let Some(control) = self
            .case
            .joint_system
            .as_ref()
            .and_then(|j| j.global_step_control.as_ref())
            .and_then(|c| c.path_error.as_ref())
        else {
            return self.trial(old, target, settings).map(|t| (t, None));
        };
        // The opt-in scope is a fixed-domain, fixed-water Boundary load path.
        // Interpolate INPUT loads and prescribed boundary factors, never solutions.
        if self.case.is_mixed_mcc()
            || target.active != old.active_elements
            || target.water_level != old.hydrostatic_water_level_m
            || target.pore_pressure_kpa != old.pore_pressure_kpa
            || target.joint_water != old.joint_water
            || target.anchor_control.is_some()
        {
            return Err(Error::Input("JOINT path error control requires a fixed domain and water field without anchor control".into()));
        }
        let start = if target.stage == old.stage_index {
            old.stage_progress
        } else {
            0.
        };
        let h = target.progress - start;
        let base = if target.stage == 0 {
            settings.initial_increment
        } else {
            self.case.stages[target.stage - 1].step(settings)
        };
        if h <= 0. || !h.is_finite() {
            return Err(Error::Input("invalid JOINT path error interval".into()));
        }
        let mut midpoint = target.clone();
        midpoint.progress = start + h * 0.5;
        if midpoint.progress == start || midpoint.progress == target.progress {
            return Err(Error::Solve(
                "JOINT path error midpoint exhausted floating point resolution".into(),
            ));
        }
        midpoint.boundary_factor =
            old.boundary_factor + (target.boundary_factor - old.boundary_factor) * 0.5;
        midpoint.pressure =
            old.surface_pressure_kpa + (target.pressure - old.surface_pressure_kpa) * 0.5;
        midpoint.load = old
            .applied_nodal_load_kn
            .iter()
            .zip(&target.load)
            .map(|(a, b)| a + (b - a) * 0.5)
            .collect();
        let old_release = if target.stage == 0 {
            &old.geostatic_release_load_kn
        } else {
            &old.excavation_release_load_kn
        };
        if old_release.len() != target.release.len() {
            return Err(Error::Input(
                "JOINT path error release dimension mismatch".into(),
            ));
        }
        midpoint.release = old_release
            .iter()
            .zip(&target.release)
            .map(|(a, b)| a + (b - a) * 0.5)
            .collect();
        let midpoint_progress = midpoint.progress;
        let coarse = self.trial(old, target.clone(), settings)?;
        let half = self.trial(old, midpoint, settings)?;
        let fine = self.trial(&half.state, target, settings)?;
        let a = old
            .joints
            .as_ref()
            .ok_or_else(|| Error::Input("missing JOINT history".into()))?;
        let b = half
            .state
            .joints
            .as_ref()
            .ok_or_else(|| Error::Input("missing midpoint JOINT history".into()))?;
        let c = coarse
            .state
            .joints
            .as_ref()
            .ok_or_else(|| Error::Input("missing coarse JOINT history".into()))?;
        let d = fine
            .state
            .joints
            .as_ref()
            .ok_or_else(|| Error::Input("missing fine JOINT history".into()))?;
        if [b.points.len(), c.points.len(), d.points.len()]
            .iter()
            .any(|n| *n != a.points.len())
        {
            return Err(Error::Input("JOINT path error point count changed".into()));
        }
        let slip_errors_m: Vec<[f64; 3]> = c
            .points
            .iter()
            .zip(&d.points)
            .map(|(c, d)| std::array::from_fn(|i| c.plastic_slip_m[i] - d.plastic_slip_m[i]))
            .collect();
        let dissipation_errors_kpa_m: Vec<f64> = c
            .points
            .iter()
            .zip(&d.points)
            .map(|(c, d)| c.plastic_dissipation_kpa_m - d.plastic_dissipation_kpa_m)
            .collect();
        // Infinity norms have a fixed GP/component traversal and no parallel reduction.
        let maximum_slip_error_m = slip_errors_m
            .iter()
            .flatten()
            .map(|x| x.abs())
            .fold(0., f64::max);
        let maximum_dissipation_error_kpa_m = dissipation_errors_kpa_m
            .iter()
            .map(|x| x.abs())
            .fold(0., f64::max);
        let slip_budget_m = (base * base) * control.slip_scale_m;
        let dissipation_budget_kpa_m = (base * base) * control.dissipation_scale_kpa_m;
        let finite = slip_errors_m
            .iter()
            .flatten()
            .chain(dissipation_errors_kpa_m.iter())
            .all(|x| x.is_finite());
        if !finite {
            return Err(Error::Solve("nonfinite JOINT path error estimate".into()));
        }
        let hidden_midpoint_switch_points = a
            .points
            .iter()
            .zip(&b.points)
            .zip(&c.points)
            .filter(|((a, b), c)| a.status == c.status && a.status != b.status)
            .count();
        let endpoint_disagreement_points = c
            .points
            .iter()
            .zip(&d.points)
            .filter(|(c, d)| c.status != d.status)
            .count();
        let audit = PathErrorAudit {
            midpoint_progress,
            base_increment: base,
            actual_increment: h,
            slip_budget_m,
            dissipation_budget_kpa_m,
            maximum_slip_error_m,
            maximum_dissipation_error_kpa_m,
            slip_errors_m,
            dissipation_errors_kpa_m,
            old_status: a.points.iter().map(|p| p.status.clone()).collect(),
            midpoint_status: b.points.iter().map(|p| p.status.clone()).collect(),
            coarse_status: c.points.iter().map(|p| p.status.clone()).collect(),
            fine_status: d.points.iter().map(|p| p.status.clone()).collect(),
            hidden_midpoint_switch_points,
            endpoint_disagreement_points,
            midpoint_relative_force_residual: half.residual,
            fine_relative_force_residual: fine.residual,
            midpoint_linear_solves: half.linear_solves,
            fine_linear_solves: fine.linear_solves,
            reject: maximum_slip_error_m > slip_budget_m
                || maximum_dissipation_error_kpa_m > dissipation_budget_kpa_m,
        };
        Ok((coarse, Some(audit)))
    }
}

/// Predict a next trial size from the O(h^2) step-doubling estimate.
/// This affects proposals only; every proposed step is solved and checked again.
pub(super) fn next_increment(a: &PathErrorAudit, minimum: f64) -> f64 {
    let ratio = (a.maximum_slip_error_m / a.slip_budget_m)
        .max(a.maximum_dissipation_error_kpa_m / a.dissipation_budget_kpa_m);
    // A fixed local error budget on a prescribed base interval: for a smooth
    // first-order path the coarse/fine difference is O(h^2).
    let factor = if ratio > 0. {
        (0.8 / ratio.sqrt()).clamp(0.5, 2.)
    } else {
        2.
    };
    (a.actual_increment * factor)
        .min(a.base_increment)
        .max(minimum)
}
