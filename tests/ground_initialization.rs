use soil_fem_3d::{
    Result,
    ground::{self, Consolidation, Ground, VoidRatio},
    mesh::{Mesh, Region},
};
use std::collections::BTreeMap;

fn fixture() -> Result<(Mesh, Ground)> {
    let mut mesh = Mesh {
        schema: "soil-fem-mesh/2".into(),
        units: "m-kN-kPa".into(),
        axis: "z_up".into(),
        crs: "LOCAL".into(),
        origin_m: [0.; 3],
        nodes: vec![],
        tetrahedra: vec![],
        tetrahedra10: vec![],
        quadrature_order: None,
        face_sets: BTreeMap::new(),
        material_ids: vec![],
        region_ids: vec![],
        excavation_sets: BTreeMap::new(),
        regions: vec![],
        node_sets: BTreeMap::new(),
    };
    for z in 0..3 {
        for y in 0..3 {
            for x in 0..3 {
                mesh.nodes
                    .push([5. * x as f64, 5. * y as f64, 5. * z as f64]);
            }
        }
    }
    for z in 0usize..2 {
        for y in 0usize..2 {
            for x in 0usize..2 {
                let id = |dx, dy, dz| (x + dx) + 3 * (y + dy) + 9 * (z + dz);
                let corners = [
                    id(0, 0, 0),
                    id(1, 0, 0),
                    id(1, 1, 0),
                    id(0, 1, 0),
                    id(0, 0, 1),
                    id(1, 0, 1),
                    id(1, 1, 1),
                    id(0, 1, 1),
                ];
                for nodes in [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
                ] {
                    let tet = nodes.map(|i| corners[i]);
                    mesh.tetrahedra.push(tet);
                    mesh.material_ids.push(z as u32 + 1);
                    mesh.region_ids.push(z as u32 + 1);
                }
            }
        }
    }
    for (i, name) in ["LowerClay", "UpperClay"].iter().enumerate() {
        mesh.regions.push(Region {
            id: i as u32 + 1,
            body: (*name).into(),
            material_id: i as u32 + 1,
            bounds_m: [[0., 0., 5. * i as f64], [10., 10., 5. * (i + 1) as f64]],
            volume_m3: 500.,
            source_sha256: "a".repeat(64),
        });
    }
    let ground = serde_json::from_str(include_str!("../examples/ground.json"))?;
    Ok((mesh, ground))
}

#[test]
fn saturated_two_layer_k0_matches_analytic_stress_and_fe_equilibrium() -> Result<()> {
    let (mesh, ground) = fixture()?;
    let state = ground::initialize(&mesh, ground)?;
    assert_eq!(state.shared_interface_triangles, 8);
    assert!(state.true_relative_free_force_residual < 1e-12);
    assert!((state.applied_force_kn[2] + 19000.).abs() < 1e-9);
    assert!((state.support_reaction_kn[2] - 19000.).abs() < 1e-9);
    for point in &state.point_states {
        let z = point.centroid_m[2];
        let sv = if z < 5. {
            90. + 20. * (5. - z)
        } else {
            18. * (10. - z)
        };
        let u = 9.80665 * (10. - z);
        assert!((point.total_stress_compression_positive_kpa[2] - sv).abs() < 1e-10);
        assert!((point.pore_pressure_kpa - u).abs() < 1e-10);
        let s = point
            .constitutive_state
            .as_mcc()
            .unwrap()
            .effective_stress_compression_positive_kpa;
        let k0 = if z < 5. { 0.6 } else { 0.5 };
        assert!((s[2] - (sv - u)).abs() < 1e-10);
        assert!((s[0] / s[2] - k0).abs() < 1e-12);
        assert!((point.mcc_pc_over_mean_effective.unwrap() - 2.).abs() < 1e-12);
        assert!(point.vertical_ocr.is_none());
        assert!(point.normalized_yield <= 0.);
    }
    Ok(())
}

#[test]
fn material_and_region_order_do_not_determine_assignment() -> Result<()> {
    let (mut mesh, mut ground) = fixture()?;
    mesh.regions.reverse();
    ground.materials.reverse();
    ground.layers.reverse();
    let result = ground::initialize(&mesh, ground)?;
    assert!((result.body_volumes_m3["LowerClay"] - 500.).abs() < 1e-9);
    assert_eq!(result.point_states[0].material_id, 1);
    Ok(())
}

#[test]
fn surcharge_and_dry_state_have_correct_applied_resultant() -> Result<()> {
    let (mesh, mut ground) = fixture()?;
    ground.initialization.water_level_m = None;
    ground.initialization.surface_surcharge_kpa = 25.;
    let state = ground::initialize(&mesh, ground)?;
    assert!((state.applied_force_kn[2] + 21500.).abs() < 1e-9);
    assert!(state.true_relative_free_force_residual < 1e-12);
    assert!(state.point_states.iter().all(|p| p.pore_pressure_kpa == 0.));
    Ok(())
}

#[test]
fn vertical_ocr_requires_pc_and_is_not_silently_an_mcc_ratio() -> Result<()> {
    let (mesh, mut ground) = fixture()?;
    ground.layers[0].consolidation = Consolidation::VerticalOcr {
        ratio: 3.,
        preconsolidation_pressure_kpa: 300.,
    };
    let result = ground::initialize(&mesh, ground)?;
    let p = &result.point_states[0];
    assert_eq!(p.vertical_ocr, Some(3.));
    assert!((p.mcc_pc_over_mean_effective.unwrap() - 3.).abs() > 1e-2);
    assert!(
        (p.past_maximum_vertical_effective_stress_kpa.unwrap_or(0.)
            - 3. * p
                .constitutive_state
                .as_mcc()
                .unwrap()
                .effective_stress_compression_positive_kpa[2])
            .abs()
            < 1e-10
    );
    let ambiguous = r#"{"definition":"vertical_ocr","ratio":3.}"#;
    assert!(serde_json::from_str::<Consolidation>(ambiguous).is_err());
    Ok(())
}

#[test]
fn rejects_yield_inconsistent_ocr_void_ratio_wrong_material_and_geometry() -> Result<()> {
    let (mesh, ground) = fixture()?;
    let mut bad = ground.clone();
    bad.layers[0].consolidation = Consolidation::MccPcOverMeanEffective { ratio: 1. };
    assert!(ground::initialize(&mesh, bad).is_err());
    let mut bad = ground.clone();
    bad.layers[0].void_ratio = VoidRatio::Explicit { value: 0.8 };
    assert!(ground::initialize(&mesh, bad).is_err());
    let mut bad = mesh.clone();
    bad.material_ids[0] = 2;
    assert!(ground::initialize(&bad, ground.clone()).is_err());
    let mut bad = mesh.clone();
    bad.regions[1].volume_m3 = 450.;
    assert!(ground::initialize(&bad, ground.clone()).is_err());
    let mut bad = ground.clone();
    bad.initialization.water_level_m = Some(11.);
    assert!(ground::initialize(&mesh, bad).is_err());
    let mut bad = ground.clone();
    bad.initialization.water_level_m = Some(7.5);
    assert!(ground::initialize(&mesh, bad).is_err());
    Ok(())
}

#[test]
fn detached_interface_and_sloping_layer_are_rejected() -> Result<()> {
    let (mesh, ground) = fixture()?;
    let mut bad = mesh.clone();
    let mut copies = BTreeMap::new();
    for i in 0..mesh.nodes.len() {
        if mesh.nodes[i][2] == 5. {
            copies.insert(i, bad.nodes.len());
            bad.nodes.push(mesh.nodes[i]);
        }
    }
    for (e, tet) in bad.tetrahedra.iter_mut().enumerate() {
        if bad.region_ids[e] == 2 {
            for node in tet {
                if let Some(&copy) = copies.get(node) {
                    *node = copy;
                }
            }
        }
    }
    let err =
        ground::initialize(&bad, ground.clone()).expect_err("non-conformal contact must fail");
    assert!(err.to_string().contains("conforming interface"));
    let mut bad = mesh;
    bad.nodes[26][2] += 1.;
    assert!(ground::initialize(&bad, ground).is_err());
    Ok(())
}
