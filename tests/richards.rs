use nalgebra::{Matrix3, Vector3};
use soil_fem_3d::{
    Result, boundary, elements,
    mesh::Mesh,
    mixed_flow,
    retention::VanGenuchten,
    richards::{self, Case},
};

fn law() -> VanGenuchten {
    VanGenuchten {
        residual_water_content: 0.08,
        saturated_water_content: 0.45,
        alpha_per_m: 1.2,
        n: 1.6,
        saturated_specific_storage_per_m: 1e-5,
    }
}

#[test]
fn finite_width_reservoir_converges_without_changing_local_mass_or_atmospheric_gates() -> Result<()>
{
    let mesh: Mesh = serde_json::from_str(include_str!("data/richards-finite-width-mesh.json"))?;
    let case: Case = serde_json::from_str(include_str!("../examples/richards-dam.json"))?;
    let report = richards::solve(&mesh, case, &vec![true; mesh.tetrahedra.len()], None)?;
    assert!(report.complete, "{:?}", report.failure_reason);
    let frame = report.frames.last().expect("steady frame");
    assert!(frame.balance.maximum_scaled_residual <= 1.);
    assert!(frame.balance.absolute_global_mass_residual_m3_s < 1e-13);
    let mut transverse_squared = 0.;
    for e in 0..mesh.tetrahedra.len() {
        let cell = mixed_flow::Cell::from_mesh(&mesh, e)?;
        let flows = frame.state.cell_outflow_m3_s[e].expect("active flux");
        for p in &cell.points {
            let qy = (0..4).map(|j| flows[j] * p.basis_per_m2[j][1]).sum::<f64>();
            transverse_squared += p.weight_m3 * qy * qy;
        }
    }
    assert!(transverse_squared.sqrt() > 1e-8);
    assert!(frame.balance.atmospheric_faces.iter().any(|f| f.seeping));
    Ok(())
}
fn relative(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}
#[test]
fn retention_matches_independent_100_digit_decimal_values_and_preserves_wet_mass_difference()
-> Result<()> {
    for line in include_str!("data/retention-decimal-reference.csv")
        .lines()
        .skip(1)
    {
        let row: Vec<f64> = line
            .split(',')
            .map(|v| v.parse().expect("decimal reference"))
            .collect();
        let s = law().evaluate(row[0])?;
        for (a, b) in [
            (s.effective_saturation, row[1]),
            (s.water_content, row[2]),
            (s.water_deficit, row[3]),
            (s.capacity_per_m, row[4]),
            (s.relative_conductivity, row[5]),
        ] {
            assert!(
                relative(a, b) < 2e-12 || (a == 0. && b == 0.),
                "head={} actual={a:e} expected={b:e}",
                row[0]
            );
        }
    }
    let a = law().evaluate(-1e-20)?;
    let b = law().evaluate(-1e-21)?;
    assert_eq!(a.water_content, law().saturated_water_content);
    assert_eq!(a.effective_saturation, 1.);
    assert!(a.water_deficit > 0. && b.storage_offset > a.storage_offset);
    let mut old = law().evaluate(-1e8)?;
    for psi in [-1e3, -100., -10., -1., -0.1, -1e-3, -1e-9, 0., 1., 100.] {
        let s = law().evaluate(psi)?;
        assert!(
            s.effective_saturation >= old.effective_saturation
                && s.relative_conductivity >= old.relative_conductivity
                && s.storage_offset >= old.storage_offset
        );
        old = s;
    }
    let step = 1e-6;
    let derivative = (law().evaluate(-0.1 + step)?.storage_offset
        - law().evaluate(-0.1 - step)?.storage_offset)
        / (2. * step);
    assert!(relative(derivative, law().evaluate(-0.1)?.capacity_per_m) < 1e-8);
    for p in [-10., -1., -0.1, -1e-4] {
        let step: f64 = -p * 1e-5;
        let derivative = (law().evaluate(p + step)?.log_relative_conductivity
            - law().evaluate(p - step)?.log_relative_conductivity)
            / (2. * step);
        assert!(
            relative(
                derivative,
                law().evaluate(p)?.d_log_relative_conductivity_per_m
            ) < 1e-7
        );
    }
    let dry = law().evaluate(-1e308)?;
    assert_eq!(dry.relative_conductivity, 0.);
    assert!(dry.log_relative_conductivity.is_finite());
    assert_eq!(law().evaluate(-1e-308)?.relative_conductivity, 1.);
    let mut bad = law();
    bad.n = 1.;
    assert!(bad.evaluate(-1.).is_err());
    assert!(law().evaluate(f64::NAN).is_err());
    Ok(())
}
fn mesh(quadratic: bool) -> Result<Mesh> {
    let name = if quadratic {
        "examples/seepage-column-tet10-mesh.json"
    } else {
        "examples/seepage-column-tet4-mesh.json"
    };
    let mut mesh: Mesh = serde_json::from_slice(&std::fs::read(name)?)?;
    let owners = boundary::surface_owners(&mesh);
    mesh.face_sets.insert(
        "exterior".into(),
        owners
            .iter()
            .filter(|(_, o)| o.len() == 1)
            .map(|(f, _)| *f)
            .collect(),
    );
    Ok(mesh)
}
fn base() -> Result<Case> {
    Ok(serde_json::from_value(serde_json::json!({
        "schema":"soil-fem-richards/1","units":"m-s",
        "materials":([1,2].iter().map(|id|serde_json::json!({"id":id,"saturated_conductivity_m_s":[[2e-6,0.2e-6,0.],[0.2e-6,1e-6,0.],[0.,0.,1.5e-6]],"retention":law()})).collect::<Vec<_>>()),
        "initial_head":{"head_m":2.},"analysis":{"kind":"steady","boundaries":[{"face_set":"exterior","condition":{"kind":"head","value":{"head_m":2.}}}]}
    }))?)
}

#[test]
fn rt0_mass_matches_exact_simplex_moments_and_piola_preserves_curved_face_fluxes() -> Result<()> {
    let linear_mesh = mesh(false)?;
    let cell = mixed_flow::Cell::from_mesh(&linear_mesh, 0)?;
    let k = Matrix3::new(2., 0.2, 0.1, 0.2, 1., 0.3, 0.1, 0.3, 3.);
    let inv = k.try_inverse().expect("SPD");
    let mass = cell.mass(&vec![inv; cell.points.len()])?;
    let v = linear_mesh
        .points(&linear_mesh.tetrahedra[0])
        .map(Vector3::from);
    let mean = v.iter().copied().sum::<Vector3<f64>>() / 4.;
    let mut moment = Matrix3::zeros();
    for i in 0..4 {
        for j in 0..4 {
            moment += v[i] * v[j].transpose() * if i == j { 0.1 } else { 0.05 };
        }
    }
    for i in 0..4 {
        for j in 0..4 {
            let exact =
                ((inv * moment).trace() - v[i].dot(&(inv * mean)) - mean.dot(&(inv * v[j]))
                    + v[i].dot(&(inv * v[j])))
                    / (9. * cell.volume_m3);
            assert!(relative(mass[(i, j)], exact) < 1e-11 || (mass[(i, j)] - exact).abs() < 1e-12);
        }
    }
    let mut curved = mesh(true)?;
    // A global quadratic shear, x'=x+0.03 z^2, is continuous and invertible.
    for p in &mut curved.nodes {
        p[0] += 0.03 * p[2] * p[2];
    }
    curved.validate()?;
    let cell = mixed_flow::Cell::from_mesh(&curved, 0)?;
    assert!(
        (cell
            .points
            .iter()
            .map(|q| q.weight_m3 * q.divergence_per_m3)
            .sum::<f64>()
            - 1.)
            .abs()
            < 2e-14
    );
    let xyz: Vec<_> = curved
        .cell_nodes(0)
        .iter()
        .map(|&n| curved.nodes[n])
        .collect();
    let rv = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].map(Vector3::from);
    for (fi, face) in elements::FACES.into_iter().enumerate() {
        let mut flux = [0.; 4];
        for weights in [
            [2. / 3., 1. / 6., 1. / 6.],
            [1. / 6., 2. / 3., 1. / 6.],
            [1. / 6., 1. / 6., 2. / 3.],
        ] {
            let r = (0..3)
                .map(|i| rv[face[i]] * weights[i])
                .sum::<Vector3<f64>>();
            let (_, j) = elements::tet_mapping(&xyz, r.into())?;
            let area =
                (j * (rv[face[1]] - rv[face[0]])).cross(&(j * (rv[face[2]] - rv[face[0]]))) / 6.;
            for i in 0..4 {
                flux[i] += (j * (r - rv[i]) * (2. / j.determinant())).dot(&area);
            }
        }
        for (i, q) in flux.into_iter().enumerate() {
            assert!(
                (q - if i == fi { 1. } else { 0. }).abs() < 3e-14,
                "face {fi} basis {i}: {q}"
            );
        }
    }
    Ok(())
}
#[test]
fn mixed_saturated_affine_tensor_flow_is_exact_on_straight_tet4_and_tet10() -> Result<()> {
    for quadratic in [false, true] {
        let mesh = mesh(quadratic)?;
        let mut c = base()?;
        let gradient = [0.2, -0.1, -0.4];
        c.initial_head.gradient = gradient;
        if let richards::Analysis::Steady { boundaries } = &mut c.analysis
            && let richards::Condition::Head { value } = &mut boundaries[0].condition
        {
            value.gradient = gradient;
        }
        let r = richards::solve(&mesh, c, &vec![true; mesh.tetrahedra.len()], None)?;
        assert!(r.complete, "{:?}", r.failure_reason);
        let final_frame = r.frames.last().expect("frame");
        let expected = -(Matrix3::new(2e-6, 0.2e-6, 0., 0.2e-6, 1e-6, 0., 0., 0., 1.5e-6)
            * Vector3::from(gradient));
        for p in &final_frame.points {
            assert!(
                (Vector3::from(p.darcy_flux_m_s) - expected).norm() < 1e-15,
                "{:?} expected {:?}",
                p.darcy_flux_m_s,
                expected
            );
        }
        assert!(final_frame.balance.maximum_absolute_cell_mass_residual_m3_s < 1e-15);
        assert!(final_frame.balance.maximum_absolute_face_flux_residual_m3_s < 1e-15);
    }
    Ok(())
}
#[test]
fn unsaturated_hydrostatic_head_has_zero_flux_and_retention_varies_with_height() -> Result<()> {
    let mut m = mesh(true)?;
    for p in &mut m.nodes {
        p[0] += 0.03 * p[2] * p[2];
    }
    let mut c = base()?;
    c.initial_head.head_m = 0.;
    if let richards::Analysis::Steady { boundaries } = &mut c.analysis
        && let richards::Condition::Head { value } = &mut boundaries[0].condition
    {
        value.head_m = 0.;
    }
    let r = richards::solve(&m, c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete, "{:?}", r.failure_reason);
    let f = r.frames.last().expect("frame");
    assert!(f.points.iter().all(|p| p.darcy_flux_m_s == [0.; 3]));
    assert!(
        f.points
            .iter()
            .any(|p| p.retention.effective_saturation < 0.8)
    );
    assert_eq!(f.balance.maximum_scaled_residual, 0.);
    Ok(())
}

fn vertical_case(bottom: f64) -> Result<Case> {
    let mut c = base()?;
    for m in &mut c.materials {
        m.saturated_conductivity_m_s = [[1e-6, 0., 0.], [0., 1e-6, 0.], [0., 0., 1e-6]];
    }
    c.analysis = serde_json::from_value(serde_json::json!({"kind":"steady","boundaries":[
        {"face_set":"bottom","condition":{"kind":"head","value":{"head_m":bottom}}},
        {"face_set":"top","condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":0.}}
    ]}))?;
    Ok(c)
}
#[test]
fn atmospheric_face_opens_for_outflow_and_closes_without_inventing_inflow() -> Result<()> {
    let m = mesh(false)?;
    for bottom in [2., 0.] {
        let c = vertical_case(bottom)?;
        let r = richards::solve(&m, c, &vec![true; m.tetrahedra.len()], None)?;
        assert!(r.complete, "bottom {bottom}: {:?}", r.failure_reason);
        let f = r.frames.last().expect("frame");
        for a in &f.balance.atmospheric_faces {
            if bottom == 2. {
                assert!(a.wet);
                assert!(a.outward_flux_m3_s > 0.);
                assert!(a.projected_pressure_head_m.abs() < 1e-9);
            } else {
                assert!(!a.wet);
                assert!(a.outward_flux_m3_s.abs() < 1e-14);
                assert!(a.projected_pressure_head_m < -0.99);
            }
        }
        let expected = if bottom == 2. { 1e-6 } else { 0. };
        let flow: f64 = f
            .balance
            .atmospheric_faces
            .iter()
            .map(|a| a.outward_flux_m3_s)
            .sum();
        assert!((flow - expected).abs() < 1e-13, "{flow}");
    }
    Ok(())
}
#[test]
fn closed_saturated_storage_gains_exactly_the_injected_volume_without_a_head_gauge() -> Result<()> {
    let m = mesh(false)?;
    let mut c = base()?;
    for mat in &mut c.materials {
        mat.volume_source_per_s = 1e-8;
    }
    c.analysis = serde_json::from_value(
        serde_json::json!({"kind":"transient","periods":[{"name":"source","end_time_s":100.,"boundaries":[]}]}),
    )?;
    c.settings.initial_time_step_s = 100.;
    c.settings.maximum_time_step_s = 100.;
    let r = richards::solve(&m, c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete, "{:?}", r.failure_reason);
    for h in r.last_committed.relative_cell_head_m.iter().flatten() {
        assert!((h - 2.1).abs() < 1e-8, "{h}");
    }
    let f = r.frames.last().expect("frame");
    assert!((f.balance.storage_rate_m3_s - 1e-8).abs() < 1e-15);
    assert!(f.balance.boundary_outflow_m3_s.abs() < 1e-15);
    assert!(f.balance.absolute_global_mass_residual_m3_s < 1e-15);
    Ok(())
}
#[test]
fn rainfall_changes_unsaturated_storage_and_failed_trials_restore_the_last_accepted_field()
-> Result<()> {
    let m = mesh(false)?;
    let mut c = vertical_case(0.)?;
    c.initial_head.head_m = 0.;
    let boundaries = serde_json::json!([
        {"face_set":"bottom","condition":{"kind":"head","value":{"head_m":0.}}},
        {"face_set":"top","condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":1e-6}}
    ]);
    c.analysis = serde_json::from_value(
        serde_json::json!({"kind":"transient","periods":[{"name":"rain","end_time_s":10000.,"boundaries":boundaries}]}),
    )?;
    c.settings.initial_time_step_s = 10000.;
    c.settings.maximum_time_step_s = 10000.;
    c.settings.maximum_iterations = 2;
    let r = richards::solve(&m, c.clone(), &vec![true; m.tetrahedra.len()], None)?;
    assert!(
        r.complete,
        "{:?}; attempts={:?}",
        r.failure_reason, r.attempts
    );
    assert!(r.frames.len() > 1);
    let f = r.frames.last().expect("frame");
    assert!(f.balance.storage_rate_m3_s > 0.);
    assert!((f.balance.precipitation_m3_s - 1e-6).abs() < 1e-18);
    let times: Vec<_> = r.frames.iter().skip(1).map(|f| f.state.time_s).collect();
    assert!(
        r.attempts.iter().any(|a| !a.accepted),
        "cutback must actually occur"
    );
    {
        c.analysis = serde_json::from_value(
            serde_json::json!({"kind":"transient","periods":times.iter().enumerate().map(|(i,t)|serde_json::json!({"name":format!("clean_{i}"),"end_time_s":t,"boundaries":boundaries})).collect::<Vec<_>>()}),
        )?;
        c.settings.maximum_iterations = 100;
        let clean = richards::solve(&m, c, &vec![true; m.tetrahedra.len()], None)?;
        assert!(clean.complete, "{:?}", clean.failure_reason);
        assert_eq!(r.last_committed, clean.last_committed);
    }
    Ok(())
}

#[test]
fn consistent_hydraulic_flux_derivative_matches_independent_directional_difference() -> Result<()> {
    let mesh = mesh(true)?;
    let cell = mixed_flow::Cell::from_mesh(&mesh, 0)?;
    let k = Matrix3::new(2e-6, 0.2e-6, 0., 0.2e-6, 1e-6, 0., 0., 0., 1.5e-6);
    let inverse = k.try_inverse().expect("SPD");
    let head = -0.2;
    let matrix = |h: f64| -> Result<nalgebra::Matrix4<f64>> {
        cell.mass(
            &cell
                .points
                .iter()
                .map(|q| Ok(inverse / law().evaluate(h - q.position_m[2])?.relative_conductivity))
                .collect::<Result<Vec<_>>>()?,
        )
    };
    let mass = matrix(head)?;
    let c = mass.try_inverse().expect("mass");
    let lambda = nalgebra::Vector4::new(-0.1, -0.3, -0.4, -0.2);
    let q = c * (nalgebra::Vector4::repeat(head) - lambda);
    let dm = cell.weighted_mass(
        &cell
            .points
            .iter()
            .map(|q| {
                let s = law().evaluate(head - q.position_m[2])?;
                Ok(-inverse * (s.d_log_relative_conductivity_per_m / s.relative_conductivity))
            })
            .collect::<Result<Vec<_>>>()?,
    )?;
    let derivative = c * nalgebra::Vector4::repeat(1.) - c * dm * q;
    let step = 1e-6;
    let flux = |h: f64| -> Result<nalgebra::Vector4<f64>> {
        Ok(matrix(h)?.try_inverse().expect("mass") * (nalgebra::Vector4::repeat(h) - lambda))
    };
    let fd = (flux(head + step)? - flux(head - step)?) / (2. * step);
    assert!((fd - derivative).norm() / fd.norm() < 1e-7);
    Ok(())
}

#[test]
fn tight_flow_regression_uses_original_residual_increments_and_conductivity_tangent() -> Result<()>
{
    let m: Mesh = serde_json::from_str(include_str!("data/richards-tight-flow-mesh.json"))?;
    let mut c: Case = serde_json::from_slice(&std::fs::read("examples/richards-rain-pause.json")?)?;
    if let richards::Analysis::Transient { periods } = &mut c.analysis {
        periods.truncate(1);
        periods[0].end_time_s = 100.;
    }
    c.settings.initial_time_step_s = 50.;
    c.settings.maximum_time_step_s = 50.;
    c.settings.absolute_flow_tolerance_m3_s = 1e-22;
    c.settings.store_integration_points = false;
    let r = richards::solve(&m, c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete, "{:?}", r.failure_reason);
    for f in r.frames.iter().skip(1) {
        assert!(f.balance.maximum_scaled_residual <= 1.);
    }
    assert_eq!(r.last_committed.time_s, 100.);
    Ok(())
}

fn checkpoint_case() -> Result<Case> {
    let mut c: Case = serde_json::from_slice(&std::fs::read("examples/richards-rain-pause.json")?)?;
    if let richards::Analysis::Transient { periods } = &mut c.analysis {
        for (i, p) in periods.iter_mut().enumerate() {
            p.end_time_s = (i + 1) as f64 * 100.;
        }
    }
    Ok(c)
}
fn normalized(r: &richards::Report) -> Result<Vec<u8>> {
    let mut r = r.clone();
    for a in &mut r.attempts {
        for i in &mut a.iterations {
            if let Some(l) = &mut i.linear {
                l.elapsed_seconds = 0.;
            }
        }
    }
    Ok(serde_json::to_vec(&r)?)
}
#[test]
fn rain_pause_restart_checkpoints_restore_all_heads_storage_faces_and_step_history() -> Result<()> {
    use soil_fem_3d::richards_checkpoint::Checkpoint;
    let m = mesh(true)?;
    let c = checkpoint_case()?;
    let active = vec![true; m.tetrahedra.len()];
    let reference = richards::solve(&m, c.clone(), &active, None)?;
    assert!(reference.complete, "{:?}", reference.failure_reason);
    let total = reference.frames.len() - 1;
    let first = reference
        .frames
        .iter()
        .position(|f| f.name == "pause")
        .expect("pause");
    let second = reference
        .frames
        .iter()
        .position(|f| f.name == "rain_again")
        .expect("rain again");
    assert!(
        reference
            .frames
            .iter()
            .any(|f| f.name == "pause" && f.balance.atmospheric_faces.iter().all(|a| !a.wet))
    );
    for count in [1, first - 1, second - 1, total - 1] {
        let partial = richards::solve_until(&m, c.clone(), &active, None, count)?;
        let checkpoint = Checkpoint::seal(partial)?;
        let decoded: Checkpoint = serde_json::from_slice(&serde_json::to_vec(&checkpoint)?)?;
        let resumed = decoded.resume(&m, c.clone(), &active, None, None)?;
        assert_eq!(normalized(&reference)?, normalized(&resumed)?);
        let mut corrupt = decoded.clone();
        corrupt.run.last_committed.cell_storage_offset_m3[0] = Some(0.);
        assert!(corrupt.resume(&m, c.clone(), &active, None, None).is_err());
        let mut altered = decoded.run.clone();
        altered.runtime.as_mut().expect("runtime").next_time_step_s *= 2.;
        let resealed = Checkpoint::seal(altered)?;
        assert!(resealed.resume(&m, c.clone(), &active, None, None).is_err());
    }
    let partial = richards::solve_until(&m, c.clone(), &active, None, 1)?;
    let saved = Checkpoint::seal(partial)?;
    let mut changed = c.clone();
    changed.materials[0].retention.alpha_per_m *= 2.;
    assert!(saved.resume(&m, changed, &active, None, None).is_err());
    let mut bad = saved.clone();
    bad.run.solver_binary_sha256 = "0".repeat(64);
    bad.payload_sha256 = soil_fem_3d::checkpoint::hash_json(&bad.run)?;
    assert!(bad.resume(&m, c, &active, None, None).is_err());
    Ok(())
}

#[test]
fn nonlinear_exhaustion_keeps_the_initial_field_and_steady_nonuniqueness_is_not_success()
-> Result<()> {
    let m = mesh(false)?;
    let active = vec![true; m.tetrahedra.len()];
    let mut c = checkpoint_case()?;
    if let richards::Analysis::Transient { periods } = &mut c.analysis {
        periods.truncate(1);
        periods[0].end_time_s = 10000.;
    }
    c.settings.initial_time_step_s = 10000.;
    c.settings.maximum_time_step_s = 10000.;
    c.settings.minimum_time_step_s = 10000.;
    c.settings.maximum_iterations = 1;
    let failed = richards::solve(&m, c, &active, None)?;
    assert!(!failed.complete && failed.failure_reason.is_some());
    assert_eq!(failed.frames.len(), 1);
    assert_eq!(failed.frames[0].state, failed.last_committed);
    assert!(
        failed
            .attempts
            .iter()
            .flat_map(|a| &a.iterations)
            .any(|i| i.line_search_fraction > 0.)
    );
    let mut c = base()?;
    c.initial_head.head_m = 0.;
    c.analysis = serde_json::from_value(
        serde_json::json!({"kind":"steady","boundaries":[{"face_set":"top","condition":{"kind":"atmospheric_no_ponding","rainfall_vertical_m_s":0.}}]}),
    )?;
    let failed = richards::solve(&m, c, &active, None)?;
    assert!(
        !failed.complete
            && failed
                .failure_reason
                .as_deref()
                .unwrap_or("")
                .contains("non-unique")
    );
    Ok(())
}

#[test]
fn hydraulic_inputs_reject_overlapping_inactive_boundaries_invalid_materials_and_insufficient_budgets()
-> Result<()> {
    let m = mesh(false)?;
    let active = vec![true; m.tetrahedra.len()];
    let original = base()?;
    let mut c = original.clone();
    c.materials[0].saturated_conductivity_m_s[0][0] = -1.;
    assert!(richards::solve(&m, c, &active, None).is_err());
    let mut c = original.clone();
    c.materials[0].retention.saturated_water_content = 0.;
    assert!(richards::solve(&m, c, &active, None).is_err());
    let mut c = original.clone();
    if let richards::Analysis::Steady { boundaries } = &mut c.analysis {
        boundaries.push(boundaries[0].clone());
    }
    assert!(richards::solve(&m, c, &active, None).is_err());
    let mut partial = active.clone();
    partial[0] = false;
    assert!(richards::solve(&m, original.clone(), &partial, None).is_err());
    assert!(richards::solve(&m, original.clone(), &active, Some(&[])).is_err());
    let mut c = vertical_case(2.)?;
    c.settings.maximum_free_faces = 1;
    let incomplete = richards::solve(&m, c, &active, None)?;
    assert!(!incomplete.complete);
    assert!(
        incomplete
            .failure_reason
            .as_deref()
            .unwrap_or("")
            .contains("budget")
    );
    assert_eq!(incomplete.frames[0].state, incomplete.last_committed);
    assert!(richards::solve_until(&m, original, &active, None, 0).is_err());
    Ok(())
}

#[test]
fn zero_storage_saturated_closed_domain_does_not_certify_an_undetermined_pressure_reference()
-> Result<()> {
    let m = mesh(false)?;
    let active = vec![true; m.tetrahedra.len()];
    let mut c = base()?;
    for material in &mut c.materials {
        material.retention.saturated_specific_storage_per_m = 0.;
    }
    c.analysis = serde_json::from_value(
        serde_json::json!({"kind":"transient","periods":[{"name":"closed","end_time_s":1.,"boundaries":[]}]}),
    )?;
    let r = richards::solve(&m, c.clone(), &active, None)?;
    assert!(
        !r.complete
            && r.failure_reason
                .as_deref()
                .unwrap_or("")
                .contains("zero-storage")
    );
    assert_eq!(r.last_committed.time_s, 0.);
    assert_eq!(r.last_committed, r.frames[0].state);
    assert!(r.attempts.iter().all(|a| !a.accepted));
    if let richards::Analysis::Transient { periods } = &mut c.analysis {
        periods[0].boundaries = serde_json::from_value(
            serde_json::json!([{"face_set":"bottom","condition":{"kind":"head","value":{"head_m":2.}}}]),
        )?;
    }
    let anchored = richards::solve(&m, c, &active, None)?;
    assert!(anchored.complete, "{:?}", anchored.failure_reason);
    assert_eq!(anchored.last_committed.time_s, 1.);
    Ok(())
}
