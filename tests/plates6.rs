use nalgebra::{DVector, Matrix3, Vector3};
use soil_fem_3d::{
    Result,
    plates::{
        Section,
        mindlin::{Kernel, shape},
    },
};
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
#[test]
fn tri6_trace_is_tet10_trace_for_all_translations_and_detects_legacy_difference() {
    let nodal = [
        [0.2, -0.1, 0.],
        [1.4, 0.2, 1e-3],
        [-0.7, 0.8, 0.00004],
        [0.9, 0.7, 0.00025],
        [-0.8, 0.6, 0.00036],
        [0.3, -0.5, 0.00001],
    ];
    let kernel = Kernel::new(points(), &section(0.1)).unwrap();
    let dofs = DVector::from_fn(36, |i, _| if i % 6 < 3 { nodal[i / 6][i % 6] } else { 0. });
    for l in [
        [0.2, 0.3, 0.5],
        [0.1, 0.8, 0.1],
        [0., 0.37, 0.63],
        [0.23, 0., 0.77],
        [0.61, 0.39, 0.],
    ] {
        let (n, _) = soil_fem_3d::elements::tet10_shape([l[1], l[2], 0.]);
        let face = [n[0], n[1], n[2], n[4], n[5], n[6]];
        for (a, plate) in kernel.translation(l, &dofs).into_iter().enumerate() {
            let solid: f64 = (0..6).map(|i| face[i] * nodal[i][a]).sum();
            assert!((solid - plate).abs() < 1e-15);
        }
    }
    let q: f64 = shape([0.5, 0.25, 0.25])
        .iter()
        .zip(nodal)
        .map(|(n, u)| n * u[2])
        .sum();
    let legacy = (nodal[3][2] + nodal[5][2]) / 2.;
    assert!((q - legacy).abs() > 1e-6);
}
#[test]
fn mindlin_six_rigid_modes_rank_and_coordinate_objectivity() -> Result<()> {
    let a = 0.43_f64;
    let r = Matrix3::new(1., 0., 0., 0., a.cos(), -a.sin(), 0., a.sin(), a.cos());
    let pts = points().map(|x| (r * Vector3::from(x) + Vector3::new(2., -3., 4.)).into());
    let k = Kernel::new(pts, &section(0.2))?;
    for mode in 0..6 {
        let rigid = DVector::from_fn(36, |i, _| {
            let p = i % 6;
            if mode < 3 {
                if p == mode { 1. } else { 0. }
            } else {
                let w = Vector3::from_fn(|a, _| if a == mode - 3 { 1. } else { 0. });
                if p < 3 {
                    w.cross(&Vector3::from(pts[i / 6]))[p]
                } else {
                    w[p - 3]
                }
            }
        });
        assert!((&k.stiffness * &rigid).norm() < 1e-12 * k.stiffness.norm() * rigid.norm());
    }
    let eig = k.stiffness.clone().symmetric_eigen().eigenvalues;
    assert_eq!(eig.iter().filter(|v| **v > 1e-10 * eig.max()).count(), 30);
    let local = Kernel::new(points(), &section(0.2))?;
    let d = DVector::from_fn(36, |i, _| ((i * 13 + 7) % 31) as f64 / 500.);
    let rd = DVector::from_fn(36, |i, _| {
        let start = (i / 3) * 3;
        (r * Vector3::new(d[start], d[start + 1], d[start + 2]))[i % 3]
    });
    let e = d.dot(&(&local.stiffness * &d));
    let re = rd.dot(&(&k.stiffness * &rd));
    assert!((e - re).abs() < 1e-11 * e.abs());
    Ok(())
}
#[test]
fn mindlin_membrane_curvature_zero_shear_and_constant_shear_patches() -> Result<()> {
    let pts = points();
    let k = Kernel::new(pts, &section(0.02))?;
    let d = DVector::from_fn(36, |i, _| {
        let [x, y, _] = pts[i / 6];
        match i % 6 {
            0 => 0.002 * x + 0.001 * y,
            1 => -0.003 * y + 0.004 * x,
            2 => -(0.005 * x * x + 0.007 * y * y) / 2. - 0.002 * x * y,
            3 => -0.007 * y - 0.002 * x,
            4 => 0.005 * x + 0.002 * y,
            _ => 0.0015,
        }
    });
    let p = k.recover("patch".into(), &d);
    for s in &p.samples {
        assert!(
            (Vector3::from(s.membrane_strain) - Vector3::new(0.002, -0.003, 0.005)).norm() < 1e-14
        );
        assert!(
            (Vector3::from(s.curvature_per_m) - Vector3::new(0.005, 0.007, 0.004)).norm() < 1e-14
        );
        assert!(s.transverse_shear_strain.iter().all(|x| x.abs() < 1e-14));
        assert!(s.drilling_mismatch_rad.abs() < 1e-14);
    }
    let energy = p.membrane_energy_kn_m
        + p.bending_energy_kn_m
        + p.transverse_shear_energy_kn_m
        + p.drilling_penalty_energy_kn_m;
    assert!((energy - 0.5 * d.dot(&(&k.stiffness * &d))).abs() < 1e-8);
    let d = DVector::from_fn(36, |i, _| {
        if i % 6 == 2 {
            0.003 * pts[i / 6][0] - 0.002 * pts[i / 6][1]
        } else {
            0.
        }
    });
    let p = k.recover("shear".into(), &d);
    for s in &p.samples {
        assert!((s.transverse_shear_strain[0] - 0.003).abs() < 1e-14);
        assert!((s.transverse_shear_strain[1] + 0.002).abs() < 1e-14);
    }
    assert!((p.transverse_shear_energy_kn_m - 0.5 * d.dot(&(&k.stiffness * &d))).abs() < 1e-8);
    Ok(())
}
#[test]
fn mindlin_consistent_traction_force_moment_and_quadratic_virtual_work_reject_equal_six()
-> Result<()> {
    let pts = points();
    let k = Kernel::new(pts, &section(0.1))?;
    let t = Vector3::new(3., -7., 2.);
    let f = k.uniform_load(t.into());
    let mut sum = Vector3::zeros();
    let mut moment = Vector3::zeros();
    for i in 0..6 {
        let fi = Vector3::new(f[6 * i], f[6 * i + 1], f[6 * i + 2]);
        sum += fi;
        moment += Vector3::from(pts[i]).cross(&fi);
    }
    let center = (Vector3::from(pts[0]) + Vector3::from(pts[1]) + Vector3::from(pts[2])) / 3.;
    assert!((sum - t * k.area_m2).norm() < 1e-12);
    assert!((moment - center.cross(&(t * k.area_m2))).norm() < 1e-12);
    assert!((0..18).all(|i| f[i] == 0.));
    // delta u_z = L0^2: integral A/6, unlike six equal nodal forces (A/4).
    let d = DVector::from_fn(36, |i, _| {
        if i % 6 == 2 {
            [1., 0., 0., 0.25, 0., 0.25][i / 6]
        } else {
            0.
        }
    });
    let work = f.dot(&d);
    assert!((work - t[2] * k.area_m2 / 6.).abs() < 1e-14);
    let equal = t[2] * k.area_m2 / 6. * d.iter().sum::<f64>();
    assert!((equal - work).abs() > 1e-3);
    Ok(())
}
#[test]
fn mindlin_rejects_curved_nonplanar_and_implicit_six_node_inputs() -> Result<()> {
    for component in 0..3 {
        let mut p = points();
        p[3][component] += 0.001;
        assert!(Kernel::new(p, &section(0.1)).is_err());
    }
    let p: soil_fem_3d::plates::Member = serde_json::from_value(
        serde_json::json!({"id":"p","nodes":["0","1","2","3","4","5"],"section_id":"s","initially_active":true}),
    )?;
    assert!(soil_fem_3d::plates::Element::new(&points(), &section(0.1), p.formulation).is_err());
    assert!(serde_json::from_value::<soil_fem_3d::plates::Member>(serde_json::json!({"id":"p","nodes":["0","1","2","3","4","5"],"section_id":"s","initially_active":true,"formulation":"mindlin_tri6","thickness_offset_m":0.1})).is_err());
    Ok(())
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
#[test]
fn mindlin_pure_bending_preserves_dkt_accuracy_across_thickness() -> Result<()> {
    let mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    for n in [1, 2, 4] {
        for h in [0.2_f64, 0.02, 0.002] {
            let model = soil_fem_3d::beams::Model::new(&mesh, strip(n, h, false)?)?;
            let ev = model.evaluate(
                &model.initial_state(),
                &vec![[0.; 3]; mesh.nodes.len()],
                &vec![true; mesh.nodes.len()],
                0.,
            )?;
            let curve = 0.002 / (2e8 * h.powi(3) / 12.);
            for (i, p) in model.positions.iter().enumerate() {
                assert!(
                    (ev.state.node_displacement_m[i][2] + curve * p[0] * p[0] / 2.).abs()
                        < 1e-7 * curve,
                    "n={n} h={h} u={:?}",
                    ev.state.node_displacement_m[i]
                );
                assert!((ev.state.node_rotation_rad[i][1] - curve * p[0]).abs() < 1e-7 * curve);
            }
            assert!(ev.report.plates.iter().all(|p|p.transverse_shear_energy_kn_m<1e-14*p.bending_energy_kn_m.max(1.)));
        }
    }
    Ok(())
}
#[test]
fn mindlin_thin_and_thick_cantilever_shear_force_converges_without_locking() -> Result<()> {
    let mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    for h in [0.8_f64, 0.02, 0.002] {
        let mut previous = 1.;
        for n in [2, 4, 8] {
            let model = soil_fem_3d::beams::Model::new(&mesh, strip(n, h, true)?)?;
            let ev = model.evaluate(
                &model.initial_state(),
                &vec![[0.; 3]; mesh.nodes.len()],
                &vec![true; mesh.nodes.len()],
                0.,
            )?;
            let exact = 0.002 * 2_f64.powi(3) / (3. * (2e8 * h.powi(3) / 12.))
                + 0.002 * 2. / ((5. / 6.) * 1e8 * h);
            let tip: Vec<_> = model
                .positions
                .iter()
                .enumerate()
                .filter(|(_, p)| p[0] == 2.)
                .map(|(i, _)| -ev.state.node_displacement_m[i][2])
                .collect();
            let mean = tip.iter().sum::<f64>() / tip.len() as f64;
            let error = (mean / exact - 1.).abs();
            assert!(
                error < 1e-7 || error < previous,
                "h={h} n={n} err={error} prev={previous}"
            );
            previous = error;
            eprintln!("cantilever h={h} n={n} relative_error={error:e}");
            if n == 8 {
                assert!(error < 1e-7, "h={h}: {error}");
            }
            assert!(ev.report.true_relative_interior_residual < 1e-10);
        }
    }
    Ok(())
}
