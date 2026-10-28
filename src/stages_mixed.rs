//! Three-field MCC Newton backend used by the existing construction controller.
//! Tet10 u, material-split continuous P1 effective mean p, independent GP volume.
//! Full system; no inversion/clipping of the material volumetric tangent.
use super::*;
use crate::{mcc, mcc_split, mixed_mcc};
#[path = "mixed_ordering.rs"]
mod ordering;

type Key = (usize, u32);

pub(super) fn maximum_material_strain(
    m: &Model<'_>,
    u: &[[f64; 3]],
    active: &[bool],
    volume: &[f64],
) -> Result<f64> {
    if u.len() != m.mesh.nodes.len()
        || active.len() != m.elements.len()
        || volume.len() != m.point_count
        || !u.iter().flatten().chain(volume).all(|v| v.is_finite())
    {
        return Err(Error::Input(
            "invalid mixed material-strain inventory or value".into(),
        ));
    }
    let mut maximum: f64 = 0.;
    for (_, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
        let ue = DVector::from_iterator(e.dofs.len(), e.dofs.iter().map(|&g| u[g / 3][g % 3]));
        for (qi, q) in e.points.iter().enumerate() {
            let strain = -&q.b * &ue;
            let correction = (volume[e.point_start + qi] - strain.iter().take(3).sum::<f64>()) / 3.;
            for i in 0..6 {
                maximum = maximum.max((strain[i] + if i < 3 { correction } else { 0. }).abs());
            }
        }
    }
    Ok(maximum)
}

fn within_initial_material_domain(
    m: &Model<'_>,
    old: &Committed,
    x: &Iterate,
    target: &Target,
) -> Result<bool> {
    if target.stage != 0 || !m.ground.is_general() {
        return Ok(true);
    }
    let cfg = m
        .ground
        .initialization
        .geostatic
        .as_ref()
        .ok_or_else(|| Error::Input("missing geostatic strain domain".into()))?;
    let history = old
        .mixed_mcc
        .as_ref()
        .ok_or_else(|| Error::Input("missing mixed strain history".into()))?;
    let cumulative: Vec<_> = history
        .cumulative_volume_strain_compression_positive
        .iter()
        .zip(&x.theta)
        .map(|(a, b)| a + b)
        .collect();
    Ok(
        maximum_material_strain(m, &x.u, &target.active, &cumulative)?
            <= cfg.maximum_initial_engineering_strain,
    )
}

pub(super) fn validate_definition(mesh: &Mesh, ground: &Ground, case: &StageCase) -> Result<()> {
    if mesh.tetrahedra10.len() != mesh.tetrahedra.len()
        || mesh.quadrature_order.unwrap_or(3) < 3
        || ground.materials.iter().any(|m| !m.model.is_mcc())
        || ground.initialization.water_level_m.is_some()
        || case.joint_system.is_some()
        || case.anchor_system.is_some()
        || !case.initially_inactive_bodies.is_empty()
        || case.stages.iter().any(|s| {
            !matches!(
                s,
                Stage::Boundary { .. }
                    | Stage::Pressure { .. }
                    | Stage::Excavate { .. }
                    | Stage::ResetDisplay { .. }
                    | Stage::InstallBeams { .. }
                    | Stage::RemoveBeams { .. }
            )
        })
    {
        return Err(Error::Input("mixed MCC requires dry Tet10/order>=3, MCC materials, and boundary/pressure/excavation/display/beam-plate installation stages; water, anchors, joints and soil material/birth transitions need explicit mixed integration".into()));
    }
    Ok(())
}

fn shape(q: &VolumePoint) -> [f64; 4] {
    [
        1. - q.natural.iter().sum::<f64>(),
        q.natural[0],
        q.natural[1],
        q.natural[2],
    ]
}

struct Layout {
    dofs: Dofs,
    pressure: BTreeMap<Key, usize>,
    theta: Vec<Option<usize>>,
    volume_rows: Vec<Vec<(usize, f64)>>,
    projected: bool,
    nt: usize,
    count: usize,
    pattern: Option<Arc<CsrPattern>>,
    ordering: Option<ordering::Ordering>,
}

impl Layout {
    fn new(m: &Model<'_>, active: &[bool], settings: &Settings, factor: f64) -> Result<Self> {
        let dofs = m.dofs(active, settings, factor)?;
        let mut theta = vec![None; m.point_count];
        let mut nt = 0;
        let mut keys = BTreeSet::new();
        for (ei, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
            for qi in 0..e.points.len() {
                theta[e.point_start + qi] = Some(nt);
                nt += 1;
            }
            for n in m.mesh.tetrahedra[ei] {
                keys.insert((n, m.mesh.material_ids[ei]));
            }
        }
        let pressure: BTreeMap<_, _> = keys.into_iter().enumerate().map(|(i, k)| (k, i)).collect();
        let projected = m.case.is_projected_mcc();
        let mut volume_rows = vec![vec![]; m.point_count];
        for (ei, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
            for (qi, q) in e.points.iter().enumerate() {
                let point = e.point_start + qi;
                volume_rows[point] = if projected {
                    m.mesh.tetrahedra[ei]
                        .iter()
                        .zip(shape(q))
                        .map(|(&n, w)| (pressure[&(n, m.mesh.material_ids[ei])], w))
                        .collect()
                } else {
                    vec![(theta[point].unwrap(), 1.)]
                };
            }
        }
        if projected {
            nt = pressure.len();
        }
        let count = nt + dofs.count + pressure.len();
        if count > settings.maximum_free_dofs || count > settings.linear_solver.dimension_limit() {
            return Err(Error::Input("mixed u/p/volume DOF budget exceeded".into()));
        }
        let pattern = if settings.linear_solver.is_dense() {
            None
        } else {
            let mut connectivity = Vec::new();
            for (ei, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
                let mut up: Vec<_> = e
                    .dofs
                    .iter()
                    .flat_map(|&g| dofs.rows[g].iter().map(|&(r, _)| nt + r))
                    .collect();
                up.extend(
                    m.mesh.tetrahedra[ei]
                        .map(|n| nt + dofs.count + pressure[&(n, m.mesh.material_ids[ei])]),
                );
                up.sort_unstable();
                up.dedup();
                if projected {
                    up.extend(
                        m.mesh.tetrahedra[ei].map(|n| pressure[&(n, m.mesh.material_ids[ei])]),
                    );
                    connectivity.push(up);
                    continue;
                }
                connectivity.push(up.clone());
                // GP volume has a diagonal block; avoid artificial GP-GP fill.
                for qi in 0..e.points.len() {
                    let t = theta[e.point_start + qi].unwrap();
                    connectivity.extend(up.iter().map(|&r| vec![t, r]));
                }
            }
            if let Some(beams) = &m.beams {
                // Condensation couples soil attachments across distinct soil
                // elements. The mixed pattern needs the same nonlocal block
                // as the displacement formulation, shifted past volume DOFs.
                let mut rows = Vec::new();
                for &node in beams.attachments.iter().flatten() {
                    for axis in 0..3 {
                        rows.extend(dofs.rows[3 * node + axis].iter().map(|&(r, _)| nt + r));
                    }
                }
                connectivity.extend(linalg::dense_block_connectivity(rows));
            }
            Some(CsrPattern::from_connectivity(count, connectivity)?)
        };
        let ordering = pattern
            .as_ref()
            .map(|p| ordering::Ordering::new(p.clone(), nt))
            .transpose()?;
        Ok(Self {
            dofs,
            pressure,
            theta,
            volume_rows,
            projected,
            nt,
            count,
            pattern,
            ordering,
        })
    }
    fn pi(&self, m: &Model<'_>, ei: usize) -> [usize; 4] {
        m.mesh.tetrahedra[ei].map(|n| self.pressure[&(n, m.mesh.material_ids[ei])])
    }
    fn matrix(&self) -> Tangent {
        self.pattern.as_ref().map_or_else(
            || Tangent::Dense(DMatrix::zeros(self.count, self.count)),
            |p| Tangent::Sparse(CsrMatrix::zeros(p.clone())),
        )
    }
}

fn pressure_projection(m: &Model<'_>, state: &Committed, layout: &Layout) -> Result<Vec<f64>> {
    let np = layout.pressure.len();
    let mut matrix = if m.case.settings.linear_solver.is_dense() {
        Tangent::Dense(DMatrix::zeros(np, np))
    } else {
        let elements = m
            .elements
            .iter()
            .enumerate()
            .filter(|(i, _)| state.active_elements[*i])
            .map(|(i, _)| layout.pi(m, i).to_vec());
        Tangent::Sparse(CsrMatrix::zeros(CsrPattern::from_connectivity(
            np, elements,
        )?))
    };
    let mut rhs = DVector::zeros(np);
    for (ei, e) in m
        .elements
        .iter()
        .enumerate()
        .filter(|(i, _)| state.active_elements[*i])
    {
        let pi = layout.pi(m, ei);
        for (qi, q) in e.points.iter().enumerate() {
            let l = shape(q);
            let p = mcc::invariants(state.material_states[e.point_start + qi].stress()).0;
            for a in 0..4 {
                rhs[pi[a]] += q.weight_m3 * l[a] * p;
                for b in 0..4 {
                    matrix.add(pi[a], pi[b], q.weight_m3 * l[a] * l[b])?;
                }
            }
        }
    }
    Ok(matrix.solve(&rhs, &m.case.settings)?.0.as_slice().to_vec())
}

fn primary(layout: &Layout, p: &[f64]) -> Vec<mixed_mcc::PressureDof> {
    layout
        .pressure
        .iter()
        .map(|(&(node_index, material_id), &i)| mixed_mcc::PressureDof {
            node_index,
            material_id,
            compression_positive_kpa: p[i],
        })
        .collect()
}

/// Work-conjugate mean stress for committed force/reaction/release evaluation.
/// Restricting the domain changes the L2 projection, never the material history.
pub(super) fn work_pressure(m: &Model<'_>, state: &Committed, active: &[bool]) -> Result<Vec<f64>> {
    let layout = Layout::new(m, active, &m.case.settings, state.boundary_factor)?;
    let p = if active == state.active_elements {
        let saved: BTreeMap<_, _> = state
            .mixed_mcc
            .as_ref()
            .ok_or_else(|| Error::Input("missing projected MCC state".into()))?
            .pressure_dofs
            .iter()
            .map(|d| ((d.node_index, d.material_id), d.compression_positive_kpa))
            .collect();
        layout
            .pressure
            .keys()
            .map(|k| {
                saved
                    .get(k)
                    .copied()
                    .ok_or_else(|| Error::Input("missing projected pressure".into()))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        let mut restricted = state.clone();
        restricted.active_elements = active.to_vec();
        pressure_projection(m, &restricted, &layout)?
    };
    let mut gp = vec![0.; m.point_count];
    for (ei, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
        let pi = layout.pi(m, ei);
        for (qi, q) in e.points.iter().enumerate() {
            gp[e.point_start + qi] = shape(q).iter().zip(pi).map(|(w, i)| w * p[i]).sum();
        }
    }
    Ok(gp)
}

fn compatibility_offsets(
    m: &Model<'_>,
    old: &Committed,
    active: &[bool],
    layout: &Layout,
) -> Result<Vec<f64>> {
    let mut result = vec![0.; layout.pressure.len()];
    if !layout.projected {
        return Ok(result);
    }
    let history = old
        .mixed_mcc
        .as_ref()
        .ok_or_else(|| Error::Input("missing projected history".into()))?;
    if active == old.active_elements {
        let saved = &history
            .projected_volume
            .as_ref()
            .ok_or_else(|| Error::Input("missing frozen compatibility offsets".into()))?
            .compatibility_offsets;
        for row in saved {
            let i = layout
                .pressure
                .get(&(row.node_index, row.material_id))
                .ok_or_else(|| Error::Input("unknown compatibility offset".into()))?;
            result[*i] = row.integral_m3;
        }
    } else {
        for (ei, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
            let pi = layout.pi(m, ei);
            let u = DVector::from_iterator(
                e.dofs.len(),
                e.dofs.iter().map(|&g| old.displacement_m[g / 3][g % 3]),
            );
            for (qi, q) in e.points.iter().enumerate() {
                let div = (&q.b * &u).iter().take(3).sum::<f64>();
                let defect =
                    div + history.cumulative_volume_strain_compression_positive[e.point_start + qi];
                for (a, w) in shape(q).into_iter().enumerate() {
                    result[pi[a]] += q.weight_m3 * w * defect;
                }
            }
        }
    }
    Ok(result)
}

fn projected_primary(
    layout: &Layout,
    theta: &[f64],
    offsets: &[f64],
    gap: f64,
) -> Option<mixed_mcc::ProjectedVolume> {
    layout.projected.then(|| mixed_mcc::ProjectedVolume {
        formulation: mixed_mcc::PROJECTED_FORMULATION.into(),
        volume_dofs: layout
            .pressure
            .iter()
            .map(|(&(node_index, material_id), &i)| mixed_mcc::VolumeDof {
                node_index,
                material_id,
                compression_positive: theta[i],
            })
            .collect(),
        compatibility_offsets: layout
            .pressure
            .iter()
            .map(
                |(&(node_index, material_id), &i)| mixed_mcc::CompatibilityOffset {
                    node_index,
                    material_id,
                    integral_m3: offsets[i],
                },
            )
            .collect(),
        maximum_gp_pressure_projection_difference_kpa: gap,
    })
}

fn metrics(
    m: &Model<'_>,
    state: &Committed,
    layout: &Layout,
    p: &[f64],
    cumulative: &[f64],
) -> (f64, f64, f64) {
    let mut volume = vec![0.; p.len()];
    let mut scale = vec![0.; p.len()];
    let mut pressure: f64 = 0.;
    let mut gap: f64 = 0.;
    let mut weak_pressure = vec![0.; p.len()];
    let mut pressure_scale = vec![0.; p.len()];
    let mut measure = 0.;
    for (ei, e) in m
        .elements
        .iter()
        .enumerate()
        .filter(|(i, _)| state.active_elements[*i])
    {
        let pi = layout.pi(m, ei);
        let u = DVector::from_iterator(
            e.dofs.len(),
            e.dofs.iter().map(|&g| state.displacement_m[g / 3][g % 3]),
        );
        for (qi, q) in e.points.iter().enumerate() {
            let point = e.point_start + qi;
            let l = shape(q);
            let epsilon = &q.b * &u;
            let div = epsilon.iter().take(3).sum::<f64>();
            let pg: f64 = (0..4).map(|a| l[a] * p[pi[a]]).sum();
            let pm = mcc::invariants(state.material_states[point].stress()).0;
            pressure = pressure.max((pg - pm).abs() / pg.abs().max(pm.abs()).max(1.));
            gap = gap.max((pg - pm).abs());
            measure += q.weight_m3;
            for a in 0..4 {
                volume[pi[a]] -= q.weight_m3 * l[a] * (div + cumulative[point]);
                scale[pi[a]] += q.weight_m3 * l[a] * (div.abs() + cumulative[point].abs());
                weak_pressure[pi[a]] += q.weight_m3 * l[a] * (pm - pg);
                pressure_scale[pi[a]] += q.weight_m3 * l[a] * pg.abs().max(pm.abs()).max(1.);
            }
        }
    }
    if let Some(projected) = &state
        .mixed_mcc
        .as_ref()
        .and_then(|s| s.projected_volume.as_ref())
    {
        for row in &projected.compatibility_offsets {
            if let Some(&i) = layout.pressure.get(&(row.node_index, row.material_id)) {
                volume[i] += row.integral_m3;
            }
        }
    }
    if layout.projected {
        pressure = norm(&weak_pressure) / norm(&pressure_scale).max(measure * 1e-12);
    }
    (
        norm(&volume) / norm(&scale).max(measure * 1e-12),
        pressure,
        gap,
    )
}

pub(super) fn initialize(m: &Model<'_>, state: &mut Committed) -> Result<()> {
    let layout = Layout::new(
        m,
        &state.active_elements,
        &m.case.settings,
        state.boundary_factor,
    )?;
    let p = pressure_projection(m, state, &layout)?;
    let cumulative = vec![0.; m.point_count];
    let (volume, pressure, gap) = metrics(m, state, &layout, &p, &cumulative);
    state.mixed_mcc = Some(mixed_mcc::State {
        tangent_method: mixed_mcc::TANGENT_METHOD.into(),
        pressure_dofs: primary(&layout, &p),
        reference_void_ratio: state
            .material_states
            .iter()
            .map(soil::State::void_ratio)
            .collect(),
        cumulative_volume_strain_compression_positive: cumulative,
        last_increment_partitions: vec![1; m.point_count],
        relative_volume_equation_residual: volume,
        relative_pressure_consistency_residual: pressure,
        projected_volume: projected_primary(
            &layout,
            &vec![0.; layout.nt],
            &vec![0.; layout.pressure.len()],
            gap,
        ),
    });
    if !m.ground.is_general() {
        validate_state(m, state)?;
    }
    Ok(())
}

pub(super) fn validate_state(m: &Model<'_>, state: &Committed) -> Result<()> {
    let mixed = state
        .mixed_mcc
        .as_ref()
        .ok_or_else(|| Error::Input("missing mixed MCC primary state".into()))?;
    if mixed.tangent_method != mixed_mcc::TANGENT_METHOD {
        return Err(Error::Input(
            "mixed MCC tangent method differs from this solver".into(),
        ));
    }
    if state.active_elements.len() != m.elements.len()
        || !state.active_elements.iter().any(|a| *a)
        || state.material_states.len() != m.point_count
        || state.displacement_m.len() != m.mesh.nodes.len()
        || mixed.cumulative_volume_strain_compression_positive.len() != m.point_count
        || mixed.reference_void_ratio.len() != m.point_count
        || mixed.last_increment_partitions.len() != m.point_count
        || mixed
            .last_increment_partitions
            .iter()
            .any(|n| *n == 0 || *n > 100_000)
        || !mixed
            .cumulative_volume_strain_compression_positive
            .iter()
            .chain(&mixed.reference_void_ratio)
            .all(|x| x.is_finite())
        || !mixed.relative_volume_equation_residual.is_finite()
        || !mixed.relative_pressure_consistency_residual.is_finite()
    {
        return Err(Error::Input(
            "invalid mixed MCC state dimensions or nonfinite values".into(),
        ));
    }
    let layout = Layout::new(
        m,
        &state.active_elements,
        &m.case.settings,
        state.boundary_factor,
    )?;
    if layout.projected != mixed.projected_volume.is_some() {
        return Err(Error::Input(
            "mixed volume formulation/state mismatch".into(),
        ));
    }
    if let Some(projected) = &mixed.projected_volume {
        if projected.formulation != mixed_mcc::PROJECTED_FORMULATION
            || projected.volume_dofs.len() != layout.nt
            || projected.compatibility_offsets.len() != layout.nt
            || !projected
                .maximum_gp_pressure_projection_difference_kpa
                .is_finite()
            || projected.maximum_gp_pressure_projection_difference_kpa < 0.
        {
            return Err(Error::Input(
                "invalid projected MCC primary inventory".into(),
            ));
        }
        let mut theta = vec![0.; layout.nt];
        let mut seen = BTreeSet::new();
        for row in &projected.volume_dofs {
            let key = (row.node_index, row.material_id);
            let i = layout
                .pressure
                .get(&key)
                .ok_or_else(|| Error::Input("unknown volume DOF".into()))?;
            if !seen.insert(key) || !row.compression_positive.is_finite() {
                return Err(Error::Input(
                    "duplicate/nonfinite projected volume DOF".into(),
                ));
            }
            theta[*i] = row.compression_positive;
        }
        seen.clear();
        for row in &projected.compatibility_offsets {
            let key = (row.node_index, row.material_id);
            if !layout.pressure.contains_key(&key)
                || !seen.insert(key)
                || !row.integral_m3.is_finite()
            {
                return Err(Error::Input(
                    "unknown/duplicate/nonfinite compatibility offset".into(),
                ));
            }
            if state.stage_index == 0 && row.integral_m3 != 0. {
                return Err(Error::Input(
                    "initial projected compatibility offset must be zero".into(),
                ));
            }
        }
        for (point, rows) in layout
            .volume_rows
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.is_empty())
        {
            let value: f64 = rows.iter().map(|&(i, w)| theta[i] * w).sum();
            if (value - mixed.cumulative_volume_strain_compression_positive[point]).abs()
                > 1e-11 * (1. + value.abs())
            {
                return Err(Error::Input(
                    "projected nodal and GP volume histories disagree".into(),
                ));
            }
        }
    }
    if mixed.pressure_dofs.len() != layout.pressure.len() {
        return Err(Error::Input(
            "mixed MCC pressure inventory differs from active mesh".into(),
        ));
    }
    let mut p = vec![0.; layout.pressure.len()];
    let mut seen = BTreeSet::new();
    for d in &mixed.pressure_dofs {
        let key = (d.node_index, d.material_id);
        let i = layout
            .pressure
            .get(&key)
            .ok_or_else(|| Error::Input("unknown mixed pressure DOF".into()))?;
        if !seen.insert(key) || !d.compression_positive_kpa.is_finite() {
            return Err(Error::Input(
                "duplicate/nonfinite mixed pressure DOF".into(),
            ));
        }
        p[*i] = d.compression_positive_kpa;
    }
    for (i, s) in state.material_states.iter().enumerate() {
        let expected = (1. + mixed.reference_void_ratio[i])
            * (-mixed.cumulative_volume_strain_compression_positive[i]).exp()
            - 1.;
        if mixed.reference_void_ratio[i] <= 0.
            || !expected.is_finite()
            || (expected - s.void_ratio()).abs() > 1e-10 * (1. + expected.abs())
        {
            return Err(Error::Input(
                "mixed cumulative volume and void-ratio history disagree".into(),
            ));
        }
    }
    let (volume, pressure, gap) = metrics(
        m,
        state,
        &layout,
        &p,
        &mixed.cumulative_volume_strain_compression_positive,
    );
    let tolerance = m.case.settings.relative_force_tolerance;
    if let Some(projected) = &mixed.projected_volume
        && (gap - projected.maximum_gp_pressure_projection_difference_kpa).abs()
            > 1e-11 * (1. + gap)
    {
        return Err(Error::Input(
            "cached projected pressure diagnostic differs".into(),
        ));
    }
    if volume > tolerance
        || pressure > tolerance
        || (volume - mixed.relative_volume_equation_residual).abs() > 1e-11
        || (pressure - mixed.relative_pressure_consistency_residual).abs() > 1e-11
    {
        return Err(Error::Input(
            "mixed pressure/volume constraints or recorded residuals disagree".into(),
        ));
    }
    Ok(())
}

pub(super) fn validate_domain_origin(
    m: &Model<'_>,
    origin: &Committed,
    current: &Committed,
) -> Result<()> {
    if !m.case.is_projected_mcc() {
        return Ok(());
    }
    let layout = Layout::new(
        m,
        &current.active_elements,
        &m.case.settings,
        current.boundary_factor,
    )?;
    let expected = compatibility_offsets(m, origin, &current.active_elements, &layout)?;
    let saved = &current
        .mixed_mcc
        .as_ref()
        .and_then(|s| s.projected_volume.as_ref())
        .ok_or_else(|| Error::Input("missing projected construction state".into()))?
        .compatibility_offsets;
    for row in saved {
        let i = layout
            .pressure
            .get(&(row.node_index, row.material_id))
            .ok_or_else(|| Error::Input("unknown construction compatibility offset".into()))?;
        if row.integral_m3 != expected[*i] {
            return Err(Error::Input(
                "frozen compatibility offset differs from construction origin".into(),
            ));
        }
    }
    Ok(())
}

struct Iterate {
    u: Vec<[f64; 3]>,
    p: Vec<f64>,
    theta: Vec<f64>,
    nodal_theta: Vec<f64>,
}

fn needed_partitions(
    m: &Model<'_>,
    old: &Committed,
    x: &Iterate,
    active: &[bool],
    settings: &Settings,
) -> Result<Vec<usize>> {
    let mut result = vec![1; m.point_count];
    let mut common = 1;
    for (_, e) in m.elements.iter().enumerate().filter(|(i, _)| active[*i]) {
        let du = DVector::from_iterator(
            e.dofs.len(),
            e.dofs
                .iter()
                .map(|&g| x.u[g / 3][g % 3] - old.displacement_m[g / 3][g % 3]),
        );
        for (qi, q) in e.points.iter().enumerate() {
            let point = e.point_start + qi;
            let de = -&q.b * &du;
            let div = de.iter().take(3).sum::<f64>();
            let magnitude = (0..6)
                .map(|i| {
                    (de[i]
                        + if i < 3 {
                            (x.theta[point] - div) / 3.
                        } else {
                            0.
                        })
                    .abs()
                })
                .sum::<f64>();
            let ratio = magnitude / settings.material_maximum_substep;
            let required = (ratio - 64. * f64::EPSILON * ratio.max(1.)).ceil().max(1.);
            if !required.is_finite() || required > 100_000. {
                return Err(Error::Material(
                    "mixed MCC fixed-partition budget exceeded".into(),
                ));
            }
            result[point] = required as usize;
            common = common.max(result[point]);
        }
    }
    // One time partition over the whole increment preserves uniform-material
    // patch fields. Per-GP counts selected by unrelated Newton overshoots
    // otherwise give different constitutive maps at identical final strains.
    for (ei, e) in m.elements.iter().enumerate() {
        if active[ei] {
            result[e.point_start..e.point_start + e.points.len()].fill(common);
        }
    }
    Ok(result)
}

struct Evaluation {
    states: Vec<soil::State>,
    beams: Option<crate::beams::Evaluation>,
    residual: DVector<f64>,
    matrix: Option<Tangent>,
    force: f64,
    volume: f64,
    pressure: f64,
    pressure_gap: f64,
}
impl Evaluation {
    fn merit(&self) -> f64 {
        self.force.max(self.volume).max(self.pressure)
    }
}

fn evaluate(
    m: &Model<'_>,
    old: &Committed,
    target: &Target,
    layout: &Layout,
    x: &Iterate,
    partitions: &[usize],
    tangent: bool,
) -> Result<Evaluation> {
    let mixed = old
        .mixed_mcc
        .as_ref()
        .ok_or_else(|| Error::Input("missing mixed MCC history".into()))?;
    let mut matrix = tangent.then(|| layout.matrix());
    let mut residual = DVector::zeros(layout.count);
    let mut internal = vec![0.; 3 * m.mesh.nodes.len()];
    let mut actual_internal = internal.clone();
    let mut states = old.material_states.clone();
    let mut scale = vec![0.; layout.pressure.len()];
    let mut volume_measure = 0.;
    let mut pressure_error: f64 = 0.;
    let mut pressure_gap: f64 = 0.;
    let mut pressure_scale = vec![0.; layout.nt];
    let offsets = compatibility_offsets(m, old, &target.active, layout)?;
    for (ei, e) in m
        .elements
        .iter()
        .enumerate()
        .filter(|(i, _)| target.active[*i])
    {
        let material = m.material_models[m.material_index(&old.element_material_ids, ei)?]
            .as_mcc()
            .ok_or_else(|| Error::Input("mixed backend requires MCC".into()))?;
        let du = DVector::from_iterator(
            e.dofs.len(),
            e.dofs
                .iter()
                .map(|&g| x.u[g / 3][g % 3] - old.displacement_m[g / 3][g % 3]),
        );
        let u = DVector::from_iterator(e.dofs.len(), e.dofs.iter().map(|&g| x.u[g / 3][g % 3]));
        let pi = layout.pi(m, ei);
        for (qi, q) in e.points.iter().enumerate() {
            let point = e.point_start + qi;
            let theta_rows = &layout.volume_rows[point];
            let l = shape(q);
            let de = -&q.b * &du;
            let div = de.iter().take(3).sum::<f64>();
            let mut dev = std::array::from_fn(|a| de[a] - if a < 3 { div / 3. } else { 0. });
            // Hydrostatic increments can leave subtraction roundoff larger
            // than the almost-zero deviatoric part. Preserve the exact trace
            // constraint instead of relaxing the split material input check.
            dev[2] = -(dev[0] + dev[1]);
            let history = old.material_states[point]
                .as_mcc()
                .ok_or_else(|| Error::Input("wrong mixed MCC history type".into()))?;
            let (update, mean, ds, t) = if tangent {
                let r = mcc_split::evaluate_implicit(
                    material,
                    history,
                    dev,
                    x.theta[point],
                    partitions[point],
                )?;
                (
                    r.update,
                    r.mean_pressure_compression_positive_kpa,
                    r.deviatoric_stress_compression_positive_kpa,
                    Some(r.tangent),
                )
            } else {
                let r = mcc_split::integrate(
                    material,
                    history,
                    dev,
                    x.theta[point],
                    partitions[point],
                )?;
                (
                    r.update,
                    r.mean_pressure_compression_positive_kpa,
                    r.deviatoric_stress_compression_positive_kpa,
                    None,
                )
            };
            let p: f64 = (0..4).map(|a| l[a] * x.p[pi[a]]).sum();
            pressure_error = pressure_error.max((mean - p).abs() / mean.abs().max(p.abs()).max(1.));
            pressure_gap = pressure_gap.max((mean - p).abs());
            let sig = SVector::<f64, 6>::from_fn(|a, _| ds[a] + if a < 3 { p } else { 0. });
            let actual = SVector::<f64, 6>::from_row_slice(
                &update.state.effective_stress_compression_positive_kpa,
            );
            let f = -q.weight_m3 * q.b.transpose() * sig;
            let actual_f = -q.weight_m3 * q.b.transpose() * actual;
            for (i, &g) in e.dofs.iter().enumerate() {
                internal[g] += f[i];
                actual_internal[g] += actual_f[i];
            }
            states[point] = soil::State::Mcc(update.state);
            for &(ti, weight) in theta_rows {
                residual[ti] += q.weight_m3 * weight * (mean - p);
                pressure_scale[ti] += q.weight_m3 * weight * mean.abs().max(p.abs()).max(1.);
            }
            let epsilon = &q.b * &u;
            let total_div = epsilon.iter().take(3).sum::<f64>();
            let total_theta =
                mixed.cumulative_volume_strain_compression_positive[point] + x.theta[point];
            volume_measure += q.weight_m3;
            for a in 0..4 {
                residual[layout.nt + layout.dofs.count + pi[a]] -=
                    q.weight_m3 * l[a] * (total_div + total_theta);
                scale[pi[a]] += q.weight_m3 * l[a] * (total_div.abs() + total_theta.abs());
            }
            if let (Some(matrix), Some(t)) = (&mut matrix, t) {
                let ke = q.weight_m3 * q.b.transpose() * t.deviatoric_by_deviatoric * &q.b;
                for &(ti, wi) in theta_rows {
                    for &(tj, wj) in theta_rows {
                        matrix.add(ti, tj, q.weight_m3 * wi * wj * t.mean_by_volume)?;
                    }
                }
                for (i, &g) in e.dofs.iter().enumerate() {
                    let bi = q.b.column(i);
                    let vol_i = bi[0] + bi[1] + bi[2];
                    let ut = -q.weight_m3
                        * (0..6)
                            .map(|a| bi[a] * t.deviatoric_by_volume[a])
                            .sum::<f64>();
                    let tu =
                        -q.weight_m3 * (0..6).map(|a| t.mean_by_deviatoric[a] * bi[a]).sum::<f64>();
                    for &(r, cr) in &layout.dofs.rows[g] {
                        let r = layout.nt + r;
                        for &(ti, wi) in theta_rows {
                            matrix.add(r, ti, cr * ut * wi)?;
                            matrix.add(ti, r, cr * tu * wi)?;
                        }
                        for (j, &h) in e.dofs.iter().enumerate() {
                            for &(c, cc) in &layout.dofs.rows[h] {
                                matrix.add(r, layout.nt + c, cr * cc * ke[(i, j)])?;
                            }
                        }
                        for a in 0..4 {
                            let pr = layout.nt + layout.dofs.count + pi[a];
                            let v = -q.weight_m3 * vol_i * l[a] * cr;
                            matrix.add(r, pr, v)?;
                            matrix.add(pr, r, v)?;
                        }
                    }
                }
                for a in 0..4 {
                    let pr = layout.nt + layout.dofs.count + pi[a];
                    for &(ti, wi) in theta_rows {
                        matrix.add(ti, pr, -q.weight_m3 * l[a] * wi)?;
                        matrix.add(pr, ti, -q.weight_m3 * l[a] * wi)?;
                    }
                }
            }
        }
    }
    let beams = if let Some(model) = &m.beams {
        let history = old
            .beams
            .as_ref()
            .ok_or_else(|| Error::Input("missing mixed beam/plate history".into()))?;
        let beam = model.evaluate(
            history,
            &x.u,
            &layout.dofs.active_nodes,
            target.boundary_factor,
        )?;
        for (i, &g) in beam.interface_dofs.iter().enumerate() {
            internal[g] += beam.interface_internal_kn[i];
            actual_internal[g] += beam.interface_internal_kn[i];
            if let Some(matrix) = &mut matrix {
                for (j, &h) in beam.interface_dofs.iter().enumerate() {
                    let value = beam.interface_tangent_kn_m[(i, j)];
                    if value != 0. {
                        for &(r, cr) in &layout.dofs.rows[g] {
                            for &(c, cc) in &layout.dofs.rows[h] {
                                matrix.add(layout.nt + r, layout.nt + c, cr * cc * value)?;
                            }
                        }
                    }
                }
            }
        }
        Some(beam)
    } else {
        None
    };
    let force: Vec<_> = internal
        .iter()
        .zip(&target.load)
        .map(|(i, f)| i - f)
        .collect();
    let actual: Vec<_> = actual_internal
        .iter()
        .zip(&target.load)
        .map(|(i, f)| i - f)
        .collect();
    let reduced = m.reduce(&force, &layout.dofs);
    for i in 0..layout.dofs.count {
        residual[layout.nt + i] = reduced[i];
    }
    let force_error = if layout.projected {
        reduced.norm()
    } else {
        reduced.norm().max(m.reduce(&actual, &layout.dofs).norm())
    } / norm(&target.load).max(1.);
    for (i, offset) in offsets.into_iter().enumerate() {
        residual[layout.nt + layout.dofs.count + i] += offset;
    }
    if layout.projected {
        pressure_error =
            residual.rows(0, layout.nt).norm() / norm(&pressure_scale).max(volume_measure * 1e-12);
    }
    let volume_error = residual
        .rows(layout.nt + layout.dofs.count, layout.pressure.len())
        .norm()
        / norm(&scale).max(volume_measure * 1e-12);
    Ok(Evaluation {
        states,
        beams,
        residual,
        matrix,
        force: force_error,
        volume: volume_error,
        pressure: pressure_error,
        pressure_gap,
    })
}

fn correction(x: &Iterate, delta: &DVector<f64>, alpha: f64, layout: &Layout) -> Iterate {
    let mut u = x.u.clone();
    for (g, row) in layout.dofs.rows.iter().enumerate() {
        for &(r, c) in row {
            u[g / 3][g % 3] += alpha * c * delta[layout.nt + r];
        }
    }
    let mut theta = x.theta.clone();
    for (i, rows) in layout.volume_rows.iter().enumerate() {
        if !rows.is_empty() {
            theta[i] += alpha * rows.iter().map(|&(r, w)| delta[r] * w).sum::<f64>();
        }
    }
    let p =
        x.p.iter()
            .enumerate()
            .map(|(i, p)| p + alpha * delta[layout.nt + layout.dofs.count + i])
            .collect();
    let nodal_theta = x
        .nodal_theta
        .iter()
        .enumerate()
        .map(|(i, t)| t + alpha * delta[i])
        .collect();
    Iterate {
        u,
        p,
        theta,
        nodal_theta,
    }
}

pub(super) fn trial(
    m: &Model<'_>,
    old: &Committed,
    target: Target,
    settings: &Settings,
) -> Result<Trial> {
    if target.pore_pressure_kpa.is_some()
        || target.water_level.is_some()
        || target.joint_water.is_some()
        || target.anchor_control.is_some()
    {
        return Err(Error::Input(
            "water/structure trial needs explicit mixed coupling".into(),
        ));
    }
    let layout = Layout::new(m, &target.active, settings, target.boundary_factor)?;
    let history = old
        .mixed_mcc
        .as_ref()
        .ok_or_else(|| Error::Input("missing mixed MCC state".into()))?;
    let saved: BTreeMap<_, _> = history
        .pressure_dofs
        .iter()
        .map(|d| ((d.node_index, d.material_id), d.compression_positive_kpa))
        .collect();
    let mut x = Iterate {
        u: old.displacement_m.clone(),
        p: vec![0.; layout.pressure.len()],
        theta: vec![0.; m.point_count],
        nodal_theta: if layout.projected {
            vec![0.; layout.nt]
        } else {
            vec![]
        },
    };
    for (key, &i) in &layout.pressure {
        x.p[i] = *saved.get(key).ok_or_else(|| {
            Error::Input("mixed pressure birth requires an explicit history policy".into())
        })?;
    }
    if let Some(t) = &layout.dofs.transform {
        let flat: Vec<_> = x.u.iter().flatten().copied().collect();
        let free = t.expand(&t.reduce(&flat));
        for (g, value) in free.iter().enumerate() {
            if layout.dofs.active_nodes[g / 3] {
                x.u[g / 3][g % 3] = t.offset[g] + value;
            }
        }
    }
    if !within_initial_material_domain(m, old, &x, &target)? {
        return Err(Error::Solve(
            "mixed geostatic material strain exceeds the configured initial small-strain domain"
                .into(),
        ));
    }
    let mut partitions = needed_partitions(m, old, &x, &target.active, settings)?;
    let mut linear_solves = Vec::new();
    for iteration in 0..=settings.maximum_newton_iterations {
        let mut eval = evaluate(m, old, &target, &layout, &x, &partitions, false)?;
        m.notify(|| {
            let mut event = m.progress(
                crate::stage_control::Kind::Newton,
                target.stage,
                target.progress,
            );
            event.newton_iteration = Some(iteration);
            event.relative_force_residual = Some(eval.force);
            event.relative_volume_residual = Some(eval.volume);
            event.relative_pressure_residual = Some(eval.pressure);
            event.maximum_material_partitions = partitions.iter().max().copied();
            event
        })?;
        if eval.merit() <= settings.relative_force_tolerance {
            let mut state = old.clone();
            state.displacement_m = x.u;
            state.material_states = eval.states;
            state.beams = eval.beams.map(|b| b.state);
            state.active_elements = target.active;
            state.applied_nodal_load_kn = target.load;
            state.boundary_factor = target.boundary_factor;
            state.surface_pressure_kpa = target.pressure;
            state.stage_index = target.stage;
            state.stage_progress = target.progress;
            state.geostatic_release_load_kn = if target.stage == 0 && m.ground.is_general() {
                target.release.clone()
            } else {
                vec![]
            };
            let transition = target.stage > 0 && m.case.stages[target.stage - 1].is_transition();
            state.excavation_release_load_kn =
                if (target.stage == 0 && m.ground.is_general()) || transition {
                    vec![0.; target.release.len()]
                } else {
                    target.release.clone()
                };
            state.transition_release_load_kn = if transition { target.release } else { vec![] };
            state.mixed_mcc = Some(mixed_mcc::State {
                tangent_method: mixed_mcc::TANGENT_METHOD.into(),
                pressure_dofs: primary(&layout, &x.p),
                reference_void_ratio: history.reference_void_ratio.clone(),
                last_increment_partitions: partitions
                    .iter()
                    .enumerate()
                    .map(|(i, n)| {
                        if layout.theta[i].is_some() {
                            *n
                        } else {
                            history.last_increment_partitions[i]
                        }
                    })
                    .collect(),
                cumulative_volume_strain_compression_positive: history
                    .cumulative_volume_strain_compression_positive
                    .iter()
                    .zip(&x.theta)
                    .map(|(a, b)| a + b)
                    .collect(),
                relative_volume_equation_residual: eval.volume,
                relative_pressure_consistency_residual: eval.pressure,
                projected_volume: if layout.projected {
                    let old_theta: BTreeMap<_, _> = history
                        .projected_volume
                        .as_ref()
                        .ok_or_else(|| Error::Input("missing nodal volume history".into()))?
                        .volume_dofs
                        .iter()
                        .map(|d| ((d.node_index, d.material_id), d.compression_positive))
                        .collect();
                    let theta = layout
                        .pressure
                        .iter()
                        .map(|(key, &i)| {
                            old_theta
                                .get(key)
                                .map(|t| t + x.nodal_theta[i])
                                .ok_or_else(|| {
                                    Error::Input(
                                        "volume birth needs explicit history policy".into(),
                                    )
                                })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    projected_primary(
                        &layout,
                        &theta,
                        &compatibility_offsets(m, old, &state.active_elements, &layout)?,
                        eval.pressure_gap,
                    )
                } else {
                    None
                },
            });
            if let Some(t) = &layout.dofs.transform {
                let error = t.maximum_constraint_error(
                    &state
                        .displacement_m
                        .iter()
                        .flatten()
                        .copied()
                        .collect::<Vec<_>>(),
                );
                if error > 1e-10 {
                    return Err(Error::Solve(
                        "mixed prescribed-displacement constraint failed".into(),
                    ));
                }
            }
            validate_state(m, &state)?;
            return Ok(Trial {
                state,
                iterations: iteration,
                residual: eval.force,
                linear_solves,
                contact_predictors: vec![],
            });
        }
        if iteration == settings.maximum_newton_iterations {
            return Err(Error::Solve(format!(
                "mixed MCC Newton limit: force={}, volume={}, pressure={}",
                eval.force, eval.volume, eval.pressure
            )));
        }
        let tangent = evaluate(m, old, &target, &layout, &x, &partitions, true)?
            .matrix
            .ok_or_else(|| Error::Solve("missing mixed tangent".into()))?;
        let rhs = -&eval.residual;
        let (delta, report) = match &layout.ordering {
            Some(ordering) => ordering.solve(&tangent, &rhs, settings)?,
            None => tangent.solve(&rhs, settings)?,
        };
        m.notify(|| {
            let mut event = m.progress(
                crate::stage_control::Kind::LinearSolved,
                target.stage,
                target.progress,
            );
            event.newton_iteration = Some(iteration);
            event.linear_solve = Some(report.clone());
            event
        })?;
        linear_solves.push(report);
        let mut alpha = 1.;
        let mut accepted = None;
        let mut repartitioned = false;
        for _ in 0..28 {
            let candidate = correction(&x, &delta, alpha, &layout);
            if !within_initial_material_domain(m, old, &candidate, &target)? {
                alpha *= 0.5;
                continue;
            }
            match needed_partitions(m, old, &candidate, &target.active, settings) {
                Ok(required) => {
                    if required.iter().zip(&partitions).any(|(a, b)| a > b) {
                        for (a, b) in partitions.iter_mut().zip(required) {
                            *a = (*a).max(b);
                        }
                        repartitioned = true;
                        break;
                    }
                    if let Ok(next) =
                        evaluate(m, old, &target, &layout, &candidate, &partitions, false)
                        && (next.merit() < eval.merit() * (1. - 1e-4 * alpha)
                            || next.merit() <= settings.relative_force_tolerance)
                    {
                        accepted = Some(candidate);
                        break;
                    }
                }
                Err(Error::Material(_)) => {}
                Err(error) => return Err(error),
            }
            alpha *= 0.5;
        }
        if repartitioned {
            continue;
        }
        x = accepted.ok_or_else(|| {
            Error::Solve(format!(
                "mixed MCC line search failed: force={}, volume={}, pressure={}",
                eval.force, eval.volume, eval.pressure
            ))
        })?;
        // Discard all trial material histories; the next evaluation restarts
        // from the unchanged committed state and the same fixed partitions.
        eval.states.clear();
    }
    unreachable!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn material_strain_includes_independent_volume_and_initial_domain_uses_cumulative_history()
    -> Result<()> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
        let ground: Ground =
            serde_json::from_str(include_str!("../examples/ground-mixed-mcc.json"))?;
        let mut case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?;
        case.stages.clear();
        let initial = solve(&mesh, ground.clone(), case.clone())?;
        let model = Model::new(&mesh, &ground, &case)?;
        let active = vec![true; mesh.tetrahedra.len()];
        let u: Vec<_> = mesh.nodes.iter().map(|p| [0., 0., -0.003 * p[2]]).collect();
        let volume = vec![0.06; model.point_count];
        assert!((maximum_material_strain(&model, &u, &active, &volume)? - 0.022).abs() < 1e-12);
        assert!((model.maximum_strain(&u, &active, false) - 0.003).abs() < 1e-12);
        let mut general: Ground =
            serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
        general
            .initialization
            .geostatic
            .as_mut()
            .unwrap()
            .maximum_initial_engineering_strain = 0.025;
        let model = Model::new(&mesh, &general, &case)?;
        let mut old = initial.last_committed;
        old.mixed_mcc
            .as_mut()
            .unwrap()
            .cumulative_volume_strain_compression_positive
            .fill(0.03);
        let mut x = Iterate {
            u,
            p: vec![],
            theta: vec![0.03; model.point_count],
            nodal_theta: vec![],
        };
        let target = Target {
            joint_water: None,
            pore_pressure_kpa: None,
            water_level: None,
            anchor_control: None,
            active,
            load: vec![],
            release: vec![],
            pressure: 100.,
            boundary_factor: 0.,
            stage: 0,
            progress: 0.5,
        };
        assert!(within_initial_material_domain(&model, &old, &x, &target)?);
        x.theta.fill(0.06); // total theta=0.09: maximum material normal strain=0.032.
        assert!(!within_initial_material_domain(&model, &old, &x, &target)?);
        let mut construction = target;
        construction.stage = 1;
        assert!(within_initial_material_domain(
            &model,
            &old,
            &x,
            &construction
        )?);
        assert!(maximum_material_strain(&model, &[], &construction.active, &volume).is_err());
        Ok(())
    }

    #[test]
    fn full_three_field_jacobian_matches_residual_directional_differences() -> Result<()> {
        check_jacobian(false, false)
    }

    #[test]
    fn projected_three_field_jacobian_matches_all_residual_blocks() -> Result<()> {
        check_jacobian(true, false)
    }

    #[test]
    fn mixed_beam_plate_jacobian_includes_nonlocal_structural_tangent() -> Result<()> {
        check_jacobian(false, true)
    }

    #[test]
    fn projected_beam_plate_jacobian_includes_nonlocal_structural_tangent() -> Result<()> {
        check_jacobian(true, true)
    }

    fn check_jacobian(projected: bool, structures: bool) -> Result<()> {
        let mesh: Mesh =
            serde_json::from_str(include_str!("../examples/hardening-column-tet10-mesh.json"))?;
        let ground: Ground =
            serde_json::from_str(include_str!("../examples/ground-mixed-mcc.json"))?;
        let mut case: StageCase =
            serde_json::from_str(include_str!("../examples/stages-mixed-mcc.json"))?;
        if projected {
            case.schema = "soil-fem-projected-mcc-stages/1".into();
        }
        if structures {
            case.beam_system = Some(serde_json::from_value(serde_json::json!({
                "schema":"soil-fem-structures/1",
                "nodes":[
                    {"id":"a","position_m":[0.,0.,1.],"match_soil_node":true},
                    {"id":"b","position_m":[1.,1.,1.],"match_soil_node":true},
                    {"id":"c","position_m":[1.,0.,1.],"match_soil_node":true}],
                "sections":[{"id":"strut","young_kpa":2e6,"poisson":0.3,"area_m2":0.02,
                    "iy_m4":0.0001,"iz_m4":0.0002,"torsion_j_m4":0.00005,"unit_weight_kn_m3":0.,
                    "formulation":{"kind":"euler_bernoulli"}}],
                "members":[{"id":"beam","nodes":["a","b"],"section_id":"strut","local_y":[0.,0.,1.],"initially_active":true}],
                "plate_sections":[{"id":"wall","young_kpa":3e6,"poisson":0.2,"thickness_m":0.05,"unit_weight_kn_m3":0.,"drilling_factor":1e-4}],
                "plates":[{"id":"plate","nodes":["a","b","c"],"section_id":"wall","initially_active":true}],
                "rotation_length_scale_m":1.,"maximum_condensed_dofs":60,
                "rotation_constraints":[
                    {"node_id":"a","direction":[1.,0.,0.],"value":0.},
                    {"node_id":"a","direction":[0.,1.,0.],"value":0.},
                    {"node_id":"a","direction":[0.,0.,1.],"value":0.}]
            }))?);
        }
        case.stages.clear();
        let initial = solve(&mesh, ground.clone(), case.clone())?;
        let old = initial.last_committed;
        let m = Model::new(&mesh, &ground, &case)?;
        let layout = Layout::new(&m, &old.active_elements, &case.settings, 0.)?;
        let mut x = Iterate {
            u: old.displacement_m.clone(),
            p: vec![0.; layout.pressure.len()],
            theta: vec![0.05; m.point_count],
            nodal_theta: if projected {
                vec![0.05; layout.nt]
            } else {
                vec![]
            },
        };
        for p in &old.mixed_mcc.as_ref().unwrap().pressure_dofs {
            x.p[layout.pressure[&(p.node_index, p.material_id)]] = p.compression_positive_kpa;
        }
        for (u, p) in x.u.iter_mut().zip(&mesh.nodes) {
            u[2] = -0.002 * p[2] * (1. - p[2]) * (1. + 0.1 * p[0]);
        }
        let target = Target {
            joint_water: None,
            pore_pressure_kpa: None,
            water_level: None,
            anchor_control: None,
            active: old.active_elements.clone(),
            load: old.applied_nodal_load_kn.clone(),
            release: old.excavation_release_load_kn.clone(),
            pressure: 100.,
            boundary_factor: 0.,
            stage: 0,
            progress: 0.,
        };
        let partitions = vec![128; m.point_count];
        let mut e = evaluate(&m, &old, &target, &layout, &x, &partitions, true)?;
        let direction = DVector::from_fn(layout.count, |i, _| {
            if i < layout.nt {
                0.001 * (1. + (i % 7) as f64 / 10.)
            } else if i < layout.nt + layout.dofs.count {
                0.001 * (i as f64).sin()
            } else {
                100. * (i as f64).cos()
            }
        });
        let predicted = e.matrix.take().unwrap().action(&direction)?;
        for h in [1e-4, 5e-5] {
            let a = correction(&x, &direction, h, &layout);
            let b = correction(&x, &direction, -h, &layout);
            let a = evaluate(&m, &old, &target, &layout, &a, &partitions, false)?;
            let b = evaluate(&m, &old, &target, &layout, &b, &partitions, false)?;
            let observed = (a.residual - b.residual) / (2. * h);
            for (start, count, tolerance) in [
                (0, layout.nt, 2e-5),
                (layout.nt, layout.dofs.count, 2e-5),
                (layout.nt + layout.dofs.count, layout.pressure.len(), 1e-7),
            ] {
                let expected = predicted.rows(start, count);
                let actual = observed.rows(start, count);
                let relative = (expected - actual).norm() / actual.norm().max(1e-12);
                assert!(
                    relative < tolerance,
                    "block start={start}, h={h}, relative={relative:e}"
                );
            }
        }
        Ok(())
    }
}
