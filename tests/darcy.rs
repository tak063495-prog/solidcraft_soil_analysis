use soil_fem_3d::{
    Result, boundary,
    darcy::{self, Case},
    mesh::Mesh,
};

fn mesh(quadratic: bool) -> Result<Mesh> {
    let path = if quadratic {
        "examples/hardening-column-tet10-mesh.json"
    } else {
        "examples/mc-column-mesh.json"
    };
    let mut m: Mesh = serde_json::from_slice(&std::fs::read(path)?)?;
    let faces = boundary::surface_owners(&m);
    for (axis, low, high) in [
        (0, "xmin", "xmax"),
        (1, "ymin", "ymax"),
        (2, "bottom", "top"),
    ] {
        for (level, name) in [(0., low), (1., high)] {
            m.face_sets.insert(
                name.into(),
                faces
                    .iter()
                    .filter(|(f, o)| {
                        o.len() == 1 && f.iter().all(|&n| (m.nodes[n][axis] - level).abs() < 1e-12)
                    })
                    .map(|(f, _)| *f)
                    .collect(),
            );
        }
    }
    Ok(m)
}
fn case(k: [[f64; 3]; 3]) -> Result<Case> {
    Ok(serde_json::from_value(serde_json::json!({
        "schema":"soil-fem-darcy/1","units":"m-s",
        "materials":[{"id":1,"hydraulic_conductivity_m_s":k},{"id":2,"hydraulic_conductivity_m_s":k}],
        "heads":[{"face_set":"bottom","head_m":2.},{"face_set":"top","head_m":1.}]
    }))?)
}
fn iso(k: f64) -> [[f64; 3]; 3] {
    [[k, 0., 0.], [0., k, 0.], [0., 0., k]]
}
#[test]
fn layered_series_flow_matches_resistance_and_nodal_balance() -> Result<()> {
    for order in [false, true] {
        let m = mesh(order)?;
        let mut c = case(iso(1e-7))?;
        c.materials[1].hydraulic_conductivity_m_s = iso(2e-6);
        let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
        assert!(r.complete, "{:?}", r.failure_reason);
        let q = 1. / (0.5 / 1e-7 + 0.5 / 2e-6);
        for (n, p) in m.nodes.iter().enumerate() {
            let h = if p[2] <= 0.5 {
                2. - q * p[2] / 1e-7
            } else {
                1. + q * (1. - p[2]) / 2e-6
            };
            assert!((r.field.head_m[n] - h).abs() < 1e-10);
        }
        assert!(
            r.points
                .iter()
                .all(|p| (p.darcy_flux_m_s[2] - q).abs() < 1e-16)
        );
        assert!(
            r.balances
                .iter()
                .all(|b| b.relative_nodal_balance_error < 1e-10)
        );
    }
    Ok(())
}
#[test]
fn quadratic_source_exact_solution_and_outflow_sign() -> Result<()> {
    let m = mesh(true)?;
    let k = 1e-6;
    let mut c = case(iso(k))?;
    c.heads[1].head_m = 2.;
    for material in &mut c.materials {
        material.volume_source_per_s = 2. * k;
    }
    let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete, "{:?}", r.failure_reason);
    for (n, p) in m.nodes.iter().enumerate() {
        assert!((r.field.head_m[n] - (2. + p[2] * (1. - p[2]))).abs() < 1e-10);
    }
    for p in &r.points {
        assert!((p.darcy_flux_m_s[2] + k * (1. - 2. * p.position_m[2])).abs() < 1e-15);
    }
    let b = &r.balances[0];
    assert!((b.source_m3_s - 2. * k).abs() < 1e-16);
    assert!((b.dirichlet_reaction_outflow_m3_s - 2. * k).abs() < 1e-16);
    assert!((b.gradient_boundary_outflow_m3_s - 2. * k).abs() < 1e-15);
    Ok(())
}
#[test]
fn rotated_anisotropic_affine_head_and_neumann_flux_on_curved_tet10() -> Result<()> {
    let mut m = mesh(true)?;
    let edge = m.tetrahedra10[0][4];
    m.nodes[edge][1] += 0.015;
    let k = [
        [2e-6, 0.3e-6, 0.2e-6],
        [0.3e-6, 1e-6, -0.1e-6],
        [0.2e-6, -0.1e-6, 3e-6],
    ];
    let gradient = [0.2, -0.1, 0.3];
    let q: [f64; 3] = std::array::from_fn(|a| -(0..3).map(|b| k[a][b] * gradient[b]).sum::<f64>());
    let mut c = case(k)?;
    c.heads.truncate(1);
    c.heads[0].face_set = "xmin".into();
    c.heads[0].head_m = 10.;
    c.heads[0].gradient = gradient;
    c.fluxes = ["xmax", "ymin", "ymax", "bottom", "top"]
        .iter()
        .map(|s| darcy::Flux::Vector {
            face_set: (*s).into(),
            flux_m_s: q,
            gradient_per_s: [[0.; 3]; 3],
            reference_point_m: [0.; 3],
        })
        .collect();
    let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete, "{:?}", r.failure_reason);
    for (n, p) in m.nodes.iter().enumerate() {
        assert!(
            (r.field.head_m[n] - (10. + (0..3).map(|a| p[a] * gradient[a]).sum::<f64>())).abs()
                < 1e-9
        );
    }
    for p in &r.points {
        for (a, expected) in q.iter().enumerate() {
            assert!((p.darcy_flux_m_s[a] - expected).abs() < 1e-14);
        }
    }
    assert!(r.balances[0].relative_nodal_balance_error < 1e-10);
    Ok(())
}
#[test]
fn saturated_solver_preserves_negative_pressure_as_incomplete() -> Result<()> {
    let m = mesh(true)?;
    let mut c = case(iso(1e-6))?;
    c.heads[0].head_m = 0.;
    c.heads[1].head_m = 0.;
    let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(!r.complete);
    assert!(r.minimum_sampled_pressure_head_m < -0.99);
    assert!(r.field.head_m.iter().all(|h| h.abs() < 1e-12));
    Ok(())
}

#[test]
fn invalid_boundaries_materials_and_inactive_faces_are_rejected() -> Result<()> {
    let mut m = mesh(true)?;
    let c = case(iso(1e-6))?;
    let active = vec![true; m.tetrahedra.len()];
    let mut no_head = c.clone();
    no_head.heads.clear();
    assert!(darcy::solve(&m, &no_head, &active, None).is_err());
    let mut bad = c.clone();
    bad.materials[0].hydraulic_conductivity_m_s[0][0] = -1e-6;
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    bad = c.clone();
    bad.materials[0].hydraulic_conductivity_m_s[0][1] = 1e-7;
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    bad = c.clone();
    bad.materials[0].hydraulic_conductivity_m_s = iso(1e308);
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    bad = c.clone();
    bad.fluxes.push(darcy::Flux::Outward {
        face_set: "top".into(),
        flux_m_s: 0.,
    });
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    m.face_sets.insert("empty".into(), vec![]);
    bad = c.clone();
    bad.heads[0].face_set = "empty".into();
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    bad = c.clone();
    bad.fluxes.push(darcy::Flux::Outward {
        face_set: "empty".into(),
        flux_m_s: 0.,
    });
    assert!(darcy::solve(&m, &bad, &active, None).is_err());
    let lower: Vec<_> = m.material_ids.iter().map(|&i| i == 1).collect();
    assert!(darcy::solve(&m, &c, &lower, None).is_err());
    assert!(darcy::solve(&m, &c, &active, Some(&vec![3; active.len()])).is_err());
    Ok(())
}
#[test]
fn malformed_head_interpolation_returns_errors_and_constant_large_head_stays_finite() -> Result<()>
{
    let field = darcy::HeadField {
        relative_head_m: vec![],
        datum_head_m: vec![],
        head_m: vec![1.],
        defined_nodes: vec![],
    };
    assert!(field.pressure_head(&[0], &[1.], 0.).is_err());
    let field = darcy::HeadField {
        relative_head_m: vec![0.],
        datum_head_m: vec![1.],
        head_m: vec![1.],
        defined_nodes: vec![true],
    };
    assert!(field.pressure_head(&[], &[], 0.).is_err());
    assert!(field.pressure_head(&[0], &[f64::NAN], 0.).is_err());
    let m = mesh(true)?;
    let mut c = case(iso(1e-6))?;
    for h in &mut c.heads {
        h.head_m = 1e308;
    }
    let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete);
    assert!(r.points.iter().all(|p| p.pressure_head_m.is_finite()));
    assert!(r.bernstein_pressure_head_lower_bound_m.is_finite());
    Ok(())
}
#[test]
fn each_disconnected_water_region_has_its_own_balance_and_needs_its_own_head() -> Result<()> {
    let mut m = mesh(true)?;
    let n = m.nodes.len();
    let count = m.tetrahedra.len();
    let added: Vec<_> = m.nodes.iter().map(|p| [p[0] + 3., p[1], p[2]]).collect();
    m.nodes.extend(added);
    m.tetrahedra
        .extend(m.tetrahedra.clone().into_iter().map(|v| v.map(|i| i + n)));
    m.tetrahedra10
        .extend(m.tetrahedra10.clone().into_iter().map(|v| v.map(|i| i + n)));
    m.material_ids.extend(m.material_ids.clone());
    m.region_ids
        .extend(m.region_ids.clone().into_iter().map(|i| i + 2));
    let mut regions = m.regions.clone();
    for r in &mut regions {
        r.id += 2;
        r.body.push_str("Second");
        for b in &mut r.bounds_m {
            b[0] += 3.;
        }
    }
    m.regions.extend(regions);
    for faces in m.face_sets.values_mut() {
        faces.extend(faces.clone().into_iter().map(|f| f.map(|i| i + n)));
    }
    let c = case(iso(1e-6))?;
    let r = darcy::solve(&m, &c, &vec![true; count * 2], None)?;
    assert!(r.complete);
    assert_eq!(r.balances.len(), 2);
    assert!(
        r.balances
            .iter()
            .all(|b| b.relative_nodal_balance_error < 1e-10)
    );
    let mut missing = m.clone();
    for faces in missing.face_sets.values_mut() {
        faces.retain(|f| f.iter().all(|&i| i < n));
    }
    assert!(darcy::solve(&missing, &c, &vec![true; count * 2], None).is_err());
    let faces = boundary::surface_owners(&m);
    m.face_sets.insert(
        "first_boundary".into(),
        faces
            .iter()
            .filter(|(f, o)| o.len() == 1 && f.iter().all(|&i| i < n))
            .map(|(f, _)| *f)
            .collect(),
    );
    m.face_sets.insert(
        "second_boundary".into(),
        faces
            .iter()
            .filter(|(f, o)| o.len() == 1 && f.iter().all(|&i| i >= n))
            .map(|(f, _)| *f)
            .collect(),
    );
    let mut zero = c.clone();
    zero.heads = serde_json::from_value(
        serde_json::json!([{"face_set":"first_boundary","head_m":2.},{"face_set":"second_boundary","head_m":3.}]),
    )?;
    let z = darcy::solve(&m, &zero, &vec![true; count * 2], None)?;
    assert!(z.complete);
    assert!(
        z.points
            .iter()
            .all(|p| p.darcy_flux_m_s.iter().all(|v| *v == 0.))
    );
    assert!(
        z.balances
            .iter()
            .all(|b| b.relative_free_head_residual == 0. && b.relative_nodal_balance_error == 0.)
    );
    zero.heads[1].head_m = 1e308;
    let huge = darcy::solve(&m, &zero, &vec![true; count * 2], None)?;
    assert!(huge.complete);
    assert!(huge.points.iter().all(|p| p.head_gradient == [0.; 3]));
    for (side, prefix) in [(false, "first"), (true, "second")] {
        for (level, suffix) in [(0., "bottom"), (1., "top")] {
            m.face_sets.insert(
                format!("{prefix}_{suffix}"),
                faces
                    .iter()
                    .filter(|(f, o)| {
                        o.len() == 1
                            && f.iter()
                                .all(|&i| (i >= n) == side && (m.nodes[i][2] - level).abs() < 1e-12)
                    })
                    .map(|(f, _)| *f)
                    .collect(),
            );
        }
    }
    for id in &mut m.material_ids[count..] {
        *id += 2;
    }
    for r in &mut m.regions[2..] {
        r.material_id += 2;
    }
    let mut contrast = c.clone();
    let mut extra = c.materials.clone();
    for material in &mut extra {
        material.id += 2;
        material.hydraulic_conductivity_m_s = iso(1e-20);
    }
    contrast.materials.extend(extra);
    contrast.heads = serde_json::from_value(
        serde_json::json!([{"face_set":"first_bottom","head_m":2.},{"face_set":"first_top","head_m":1.},{"face_set":"second_bottom","head_m":5.},{"face_set":"second_top","head_m":4.}]),
    )?;
    let result = darcy::solve(&m, &contrast, &vec![true; count * 2], None)?;
    assert!(result.complete);
    for (i, p) in m.nodes.iter().enumerate() {
        assert!((result.field.head_m[i] - ((if i < n { 2. } else { 5. }) - p[2])).abs() < 1e-10);
    }
    assert!(
        result
            .balances
            .iter()
            .all(|b| b.relative_free_head_residual < 1e-11)
    );
    Ok(())
}

#[test]
fn dense_and_sparse_backends_match_the_original_source_problem() -> Result<()> {
    let m = mesh(true)?;
    let mut c = case(iso(1e-6))?;
    c.heads[1].head_m = 2.;
    for material in &mut c.materials {
        material.volume_source_per_s = 2e-6;
    }
    let active = vec![true; m.tetrahedra.len()];
    let sparse = darcy::solve(&m, &c, &active, None)?;
    c.settings.linear_solver = soil_fem_3d::linalg::LinearSolver::DenseLu;
    c.settings.maximum_free_heads = 2000;
    let dense = darcy::solve(&m, &c, &active, None)?;
    assert!(sparse.complete && dense.complete);
    for (a, b) in sparse.field.head_m.iter().zip(&dense.field.head_m) {
        assert!((a - b).abs() < 1e-12);
    }
    assert!(dense.linear.as_ref().unwrap().true_relative_residual < 1e-11);
    Ok(())
}

#[test]
fn original_rows_are_refined_after_prescribed_head_elimination_on_a_cad_mesh() -> Result<()> {
    let m: Mesh = serde_json::from_str(include_str!("../examples/darcy-residual-mesh.json"))?;
    let c: Case = serde_json::from_str(include_str!("../examples/darcy-residual-flow.json"))?;
    let r = darcy::solve(&m, &c, &vec![true; m.tetrahedra.len()], None)?;
    assert!(r.complete);
    assert!(!r.refinement_solves.is_empty());
    assert!(r.relative_free_head_residual <= c.settings.relative_linear_tolerance);
    let flux = 0.5 / (4. / 1e-6 + 4. / 3e-6);
    for (n, p) in m.nodes.iter().enumerate() {
        let z = p[2];
        let expected = if z <= 4. {
            8.5 - flux * z / 1e-6
        } else {
            8.125 - flux * (z - 4.) / 3e-6
        };
        assert!((r.field.head_m[n] - expected).abs() < 1e-9);
    }
    Ok(())
}

#[test]
fn balanced_neumann_and_volume_exchange_use_absolute_flow_for_water_balance() -> Result<()> {
    let m = mesh(true)?;
    let active = vec![true; m.tetrahedra.len()];
    let mut c = case(iso(1e-6))?;
    c.heads.truncate(1);
    c.heads[0].gradient = [0.3, 0., 0.];
    c.fluxes = vec![
        darcy::Flux::Outward {
            face_set: "xmin".into(),
            flux_m_s: 3e-7,
        },
        darcy::Flux::Outward {
            face_set: "xmax".into(),
            flux_m_s: -3e-7,
        },
    ];
    let circulation = darcy::solve(&m, &c, &active, None)?;
    assert!(circulation.complete);
    let balance = &circulation.balances[0];
    assert!((balance.absolute_neumann_flow_m3_s - 6e-7).abs() < 1e-19);
    assert!(balance.relative_nodal_balance_error < 1e-10);
    for (n, p) in m.nodes.iter().enumerate() {
        assert!((circulation.field.head_m[n] - (2. + 0.3 * p[0])).abs() < 1e-10);
    }
    c.heads[0].gradient = [0.; 3];
    c.fluxes.clear();
    c.materials[0].volume_source_per_s = 1e-6;
    c.materials[1].volume_source_per_s = -1e-6;
    let sources = darcy::solve(&m, &c, &active, None)?;
    assert!(sources.complete);
    let balance = &sources.balances[0];
    assert!((balance.absolute_source_m3_s - 1e-6).abs() < 1e-19);
    assert!(balance.relative_nodal_balance_error < 1e-10);
    for (n, p) in m.nodes.iter().enumerate() {
        let z = p[2];
        let expected = if z <= 0.5 {
            2. - z * z / 2.
        } else {
            2.25 - z + z * z / 2.
        };
        assert!((sources.field.head_m[n] - expected).abs() < 1e-10);
    }
    let mut shifted = c.clone();
    shifted
        .materials
        .iter_mut()
        .for_each(|m| m.volume_source_per_s = 0.);
    shifted.heads[0].head_m = 1e12;
    shifted.heads[0].gradient = [0.3, 0., 0.];
    shifted.fluxes = vec![
        darcy::Flux::Outward {
            face_set: "xmin".into(),
            flux_m_s: 3e-7,
        },
        darcy::Flux::Outward {
            face_set: "xmax".into(),
            flux_m_s: -3e-7,
        },
    ];
    let large = darcy::solve(&m, &shifted, &active, None)?;
    assert!(large.complete);
    for p in &large.points {
        assert!((p.darcy_flux_m_s[0] + 3e-7).abs() < 1e-16);
    }
    Ok(())
}
