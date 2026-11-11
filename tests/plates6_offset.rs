use nalgebra::{DVector, Matrix3, Vector3};
use soil_fem_3d::{
    Result,
    plates::{Section, mindlin, offset::Kernel},
};
#[test]
fn offset_checkpoint_replays_and_refuses_rehashed_type_offset_axes_and_reference() -> Result<()> {
    use soil_fem_3d::{
        checkpoint::{Checkpoint, hash_json},
        stages,
    };
    let mesh = serde_json::from_str(include_str!("../examples/plate6-offset-mesh.json"))?;
    let ground = serde_json::from_str(include_str!("../examples/plate6-offset-ground.json"))?;
    let case = serde_json::from_str::<stages::StageCase>(include_str!(
        "../examples/plate6-offset-stages.json"
    ))?;
    let run = stages::solve_until(&mesh, ground, case, 17)?;
    let cp = Checkpoint::seal(&mesh, run)?;
    cp.verify(&mesh, &cp.run.ground, &cp.run.case)?;
    for mode in 0..4 {
        let mut bad: Checkpoint = serde_json::from_value(serde_json::to_value(&cp)?)?;
        match mode {
            0 => {
                bad.run.case.beam_system.as_mut().unwrap().plates[0].formulation =
                    soil_fem_3d::plates::Formulation::MindlinTri6
            }
            1 => bad.run.case.beam_system.as_mut().unwrap().plates[0].attachment_offset_m = -0.1,
            2 => {
                bad.run
                    .frames
                    .iter_mut()
                    .find(|f| !f.beam_results.as_ref().unwrap().plates.is_empty())
                    .unwrap()
                    .beam_results
                    .as_mut()
                    .unwrap()
                    .plates[0]
                    .local_axes[0][0] += 0.01
            }
            _ => {
                bad.run
                    .last_committed
                    .beams
                    .as_mut()
                    .unwrap()
                    .plate_reference_global_kinematics[0][0] += 0.001
            }
        }
        bad.run_sha256 = hash_json(&bad.run)?;
        assert_eq!(bad.run_sha256, hash_json(&bad.run)?);
        let error = bad
            .verify(&mesh, &cp.run.ground, &cp.run.case)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("accepted-path replay mismatch"),
            "{mode}: {error}"
        );
    }
    Ok(())
}
#[test]
fn offset_zero_is_bit_identical_to_planar_kernel() -> Result<()> {
    let base = mindlin::Kernel::new(points(), &section(0.2))?;
    let new = Kernel::new(points(), &section(0.2), 0.)?;
    let d = DVector::from_fn(36, |i, _| ((i * 17) % 23) as f64 * 1e-7);
    assert_eq!(base.stiffness, new.stiffness);
    assert_eq!(
        base.uniform_load([1., 2., 3.]),
        new.uniform_load([1., 2., 3.])
    );
    assert_eq!(
        serde_json::to_value(base.recover("p".into(), &d))?,
        serde_json::to_value(new.recover("p".into(), &d))?
    );
    Ok(())
}
#[test]
fn offset_global_rotation_and_normal_sign_are_objective() -> Result<()> {
    let angle = 0.43_f64;
    let r = Matrix3::new(
        1.,
        0.,
        0.,
        0.,
        angle.cos(),
        -angle.sin(),
        0.,
        angle.sin(),
        angle.cos(),
    );
    let pts = points().map(|p| (r * Vector3::from(p) + Vector3::new(2., -3., 4.)).into());
    let k = Kernel::new(points(), &section(0.2), 0.1)?;
    let rotated = Kernel::new(pts, &section(0.2), 0.1)?;
    let d = DVector::from_fn(36, |i, _| ((i * 13 + 7) % 31) as f64 / 500.);
    let rd = DVector::from_fn(36, |i, _| {
        let s = (i / 3) * 3;
        (r * Vector3::new(d[s], d[s + 1], d[s + 2]))[i % 3]
    });
    let energy = d.dot(&(&k.stiffness * &d));
    assert!((energy - rd.dot(&(&rotated.stiffness * &rd))).abs() < 1e-11 * energy.abs());
    let order = [0, 2, 1, 5, 4, 3];
    let reverse = Kernel::new(order.map(|i| points()[i]), &section(0.2), -0.1)?;
    let reverse_d = DVector::from_fn(36, |i, _| d[6 * order[i / 6] + i % 6]);
    assert!(
        (energy - reverse_d.dot(&(&reverse.stiffness * &reverse_d))).abs() < 1e-11 * energy.abs()
    );
    Ok(())
}
#[test]
fn offset_rigid_modes_rank_patch_energy_and_trace() -> Result<()> {
    let pts = points();
    for e in [-0.1, 0., 0.1] {
        let k = Kernel::new(pts, &section(0.2), e)?;
        for mode in 0..6 {
            let w = Vector3::from_fn(|a, _| if mode == a + 3 { 1. } else { 0. });
            let d = DVector::from_fn(36, |i, _| {
                if i % 6 < 3 {
                    w.cross(&Vector3::from(pts[i / 6]))[i % 6] + if mode == i % 6 { 1. } else { 0. }
                } else {
                    w[i % 6 - 3]
                }
            });
            assert!((&k.stiffness * &d).norm() < 1e-12 * k.stiffness.norm() * d.norm());
            let r = k.recover("rigid".into(), &d);
            for s in r.samples {
                assert!(
                    s.membrane_strain
                        .iter()
                        .chain(&s.curvature_per_m)
                        .chain(&s.transverse_shear_strain)
                        .all(|v| v.abs() < 1e-13)
                );
            }
        }
        let eig = k.stiffness.clone().symmetric_eigen().eigenvalues;
        assert_eq!(eig.iter().filter(|v| **v > 1e-10 * eig.max()).count(), 30);
        // Prescribe surface strain and curvature; neutral strain is eps_s-e*kappa.
        let d = DVector::from_fn(36, |i, _| {
            let [x, y, _] = pts[i / 6];
            match i % 6 {
                0 => 0.002 * x,
                1 => -0.003 * y,
                2 => -(0.005 * x * x + 0.007 * y * y) / 2.,
                3 => -0.007 * y,
                4 => 0.005 * x,
                _ => 0.,
            }
        });
        let r = k.recover("patch".into(), &d);
        for s in &r.samples {
            assert!(
                (Vector3::from(s.membrane_strain)
                    - Vector3::new(0.002 - e * 0.005, -0.003 - e * 0.007, 0.))
                .norm()
                    < 1e-14
            );
            assert!(
                (Vector3::from(s.curvature_per_m) - Vector3::new(0.005, 0.007, 0.)).norm() < 1e-14
            );
            assert!(s.transverse_shear_strain.iter().all(|v| v.abs() < 1e-14));
        }
        let energy = r.membrane_energy_kn_m
            + r.bending_energy_kn_m
            + r.transverse_shear_energy_kn_m
            + r.drilling_penalty_energy_kn_m;
        assert!((energy - 0.5 * d.dot(&(&k.stiffness * &d))).abs() < 1e-8);
        let arbitrary = DVector::from_fn(36, |i, _| ((i * 17) % 23) as f64 * 1e-5);
        for l in [
            [0.2, 0.3, 0.5],
            [0., 0.37, 0.63],
            [0.23, 0., 0.77],
            [0.61, 0.39, 0.],
        ] {
            let (n, _) = soil_fem_3d::elements::tet10_shape([l[1], l[2], 0.]);
            let face = [n[0], n[1], n[2], n[4], n[5], n[6]];
            let u = k.attachment_translation(l, &arbitrary);
            for a in 0..3 {
                let soil: f64 = (0..6).map(|i| face[i] * arbitrary[6 * i + a]).sum();
                assert!((u[a] - soil).abs() < 1e-10);
            }
        }
    }
    Ok(())
}
#[test]
fn offset_load_virtual_work_and_neutral_plane_moment() -> Result<()> {
    for e in [-0.1, 0.1] {
        let k = Kernel::new(points(), &section(0.2), e)?;
        let t = Vector3::new(3., -7., 2.);
        let f = k.uniform_load(t.into());
        let mut sum = Vector3::zeros();
        let mut moment = Vector3::zeros();
        for (i, p) in points().iter().enumerate() {
            let fi = f.rows(6 * i, 3).into_owned();
            sum += &fi;
            moment += Vector3::from(*p).cross(&fi) + f.rows(6 * i + 3, 3);
        }
        let center =
            (Vector3::from(points()[0]) + Vector3::from(points()[1]) + Vector3::from(points()[2]))
                / 3.
                - Vector3::new(0., 0., e);
        assert!((sum - t * k.neutral.area_m2).norm() < 1e-12);
        assert!((moment - center.cross(&(t * k.neutral.area_m2))).norm() < 1e-12);
        let d = DVector::from_fn(36, |i, _| ((i * 7) % 23) as f64 / 500.);
        assert!(
            (f.dot(&d)
                - k.neutral
                    .uniform_load(t.into())
                    .dot(&k.neutral_kinematics(&d)))
            .abs()
                < 1e-12
        );
    }
    Ok(())
}
#[test]
fn offset_cantilever_and_simply_supported_eccentric_axial_load_analytic() -> Result<()> {
    use serde_json::json;
    let mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    for simple in [false, true] {
        for n in [1, 2, 4] {
            for h in [0.2_f64, 0.02] {
                for sign in [-1., 0., 1.] {
                    let e = sign * h / 2.;
                    let force = 0.002;
                    let curve = e * force / (2e8 * h.powi(3) / 12.);
                    let strain = force / (2e8 * h) + e * curve;
                    let mut cfg = strip(n, h, false)?;
                    for p in &mut cfg.plates {
                        p.formulation = soil_fem_3d::plates::Formulation::MindlinTri6Offset;
                        p.attachment_offset_m = e;
                    }
                    cfg.nodal_loads.clear();
                    for node in &cfg.nodes {
                        let [x, y, _] = node.position_m;
                        if x == 2. {
                            let w = if y == 0.5 { 4. / 6. } else { 1. / 6. };
                            cfg.nodal_loads.push(serde_json::from_value(
                                json!({"node_id":node.id,"constant_force_kn":[force*w,0.,0.]}),
                            )?);
                        }
                    }
                    if simple {
                        cfg.rotation_constraints.clear();
                        for node in &cfg.nodes {
                            if node.position_m[0] == 2. {
                                cfg.translation_constraints.push(serde_json::from_value(
                                    json!({"node_id":node.id,"direction":[0.,0.,1.],"value":0.}),
                                )?);
                            }
                        }
                    }
                    let model = soil_fem_3d::beams::Model::new(&mesh, cfg)?;
                    let ev = model.evaluate(
                        &model.initial_state(),
                        &vec![[0.; 3]; mesh.nodes.len()],
                        &vec![true; mesh.nodes.len()],
                        0.,
                    )?;
                    assert!(ev.report.true_relative_interior_residual < 1e-10);
                    for (i, p) in model.positions.iter().enumerate() {
                        let x = p[0];
                        let expected_w = if simple {
                            curve * x * (2. - x) / 2.
                        } else {
                            -curve * x * x / 2.
                        };
                        let expected_rotation = curve * (x - if simple { 1. } else { 0. });
                        assert!(
                            (ev.state.node_displacement_m[i][0] - strain * x).abs()
                                < 1e-7 * strain.abs()
                        );
                        assert!(
                            (ev.state.node_displacement_m[i][2] - expected_w).abs()
                                < 1e-7 * curve.abs().max(strain.abs())
                        );
                        assert!(
                            (ev.state.node_rotation_rad[i][1] - expected_rotation).abs() < 1e-9
                        );
                    }
                }
            }
        }
    }
    Ok(())
}
#[test]
fn offset_rejects_outside_scope_and_requires_explicit_selection() -> Result<()> {
    use soil_fem_3d::plates::{Element, Formulation};
    for e in [f64::NAN, f64::INFINITY, 0.10001, -0.10001] {
        assert!(Kernel::new(points(), &section(0.2), e).is_err());
    }
    for a in 0..3 {
        let mut p = points();
        p[3][a] += 0.001;
        assert!(Kernel::new(p, &section(0.2), 0.01).is_err());
    }
    assert!(
        Element::new_with_offset(&points(), &section(0.2), Formulation::MindlinTri6, 0.01).is_err()
    );
    assert!(
        Element::new_with_offset(&points()[..3], &section(0.2), Formulation::CstDkt, 0.01).is_err()
    );
    let mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let mut cfg = strip(1, 0.2, false)?;
    for p in &mut cfg.plates {
        p.formulation = Formulation::MindlinTri6Offset;
        p.attachment_offset_m = 0.1;
    }
    cfg.plates[1].attachment_offset_m = -0.1;
    assert!(soil_fem_3d::beams::Model::new(&mesh, cfg).is_err());
    Ok(())
}
fn section(h: f64) -> Section {
    Section {
        id: "s".into(),
        young_kpa: 2e8,
        poisson: 0.,
        thickness_m: h,
        unit_weight_kn_m3: 0.,
        drilling_factor: 1e-4,
    }
}
fn points() -> [[f64; 3]; 6] {
    [
        [0., 0., 0.],
        [2., 0., 0.],
        [0.4, 1.3, 0.],
        [1., 0., 0.],
        [1.2, 0.65, 0.],
        [0.2, 0.65, 0.],
    ]
}
fn strip(n: usize, h: f64, shear_load: bool) -> Result<soil_fem_3d::beams::Configuration> {
    use serde_json::json;
    use std::collections::BTreeMap;
    let mut pts = Vec::<[f64; 3]>::new();
    let mut ids = BTreeMap::new();
    let mut plates = vec![];
    for i in 0..n {
        for corners in [
            [(2 * i, 0), (2 * i + 2, 0), (2 * i + 2, 2)],
            [(2 * i, 0), (2 * i + 2, 2), (2 * i, 2)],
        ] {
            let mut cell = vec![];
            for (x, y) in corners.into_iter().chain(
                [
                    (corners[0].0 + corners[1].0) / 2,
                    (corners[1].0 + corners[2].0) / 2,
                    (corners[2].0 + corners[0].0) / 2,
                ]
                .into_iter()
                .zip([
                    (corners[0].1 + corners[1].1) / 2,
                    (corners[1].1 + corners[2].1) / 2,
                    (corners[2].1 + corners[0].1) / 2,
                ]),
            ) {
                let next = pts.len();
                let node = *ids.entry((x, y)).or_insert_with(|| {
                    pts.push([x as f64 / n as f64, y as f64 / 2., 0.]);
                    next
                });
                cell.push(format!("n{node}"));
            }
            plates.push(json!({"id":format!("p{}",plates.len()),"nodes":cell,"section_id":"s","initially_active":true,"formulation":"mindlin_tri6"}));
        }
    }
    let mut translations = vec![];
    let mut rotations = vec![];
    let mut loads = vec![];
    for (i, p) in pts.iter().enumerate() {
        if p[0] == 0. {
            for a in 0..3 {
                let direction: [f64; 3] = std::array::from_fn(|b| if a == b { 1. } else { 0. });
                let c = json!({"node_id":format!("n{i}"),"direction":direction,"value":0.});
                translations.push(c.clone());
                rotations.push(c);
            }
        }
        if p[0] == 2. {
            let weight = if p[1] == 0.5 { 4. / 6. } else { 1. / 6. };
            loads.push(if shear_load {
                json!({"node_id":format!("n{i}"),"constant_force_kn":[0.,0.,-0.002*weight]})
            } else {
                json!({"node_id":format!("n{i}"),"constant_moment_kn_m":[0.,0.002*weight,0.]})
            });
        }
    }
    Ok(serde_json::from_value(
        json!({"schema":"soil-fem-structures/1","nodes":pts.iter().enumerate().map(|(i,p)|json!({"id":format!("n{i}"),"position_m":p,"soil_node_set":null})).collect::<Vec<_>>(),"sections":[],"members":[],"plate_sections":[section(h)],"plates":plates,"rotation_length_scale_m":2.,"maximum_condensed_dofs":2000,"translation_constraints":translations,"rotation_constraints":rotations,"nodal_loads":loads}),
    )?)
}
