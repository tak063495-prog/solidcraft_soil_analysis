use soil_fem_3d::{
    Result,
    checkpoint::Checkpoint,
    geostatic::{self, MccAffineField, Seed},
    ground::{Ground, VoidRatio},
    mesh::Mesh,
    stages::{self, StageCase},
};

fn fixture(quadratic: bool) -> Result<(Mesh, Ground, StageCase)> {
    let mut mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    mesh.face_sets.insert(
        "top".into(),
        soil_fem_3d::boundary::surface_owners(&mesh)
            .into_iter()
            .filter(|(f, o)| o.len() == 1 && f.iter().all(|&n| mesh.nodes[n][2] == 1.))
            .map(|(f, _)| f)
            .collect(),
    );
    let ground =
        serde_json::from_str(include_str!("../examples/ground-geostatic-mcc-column.json"))?;
    let case = serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    Ok((
        if quadratic {
            mesh.elevate_straight()?
        } else {
            mesh
        },
        ground,
        case,
    ))
}

fn fields(ground: &Ground) -> Vec<MccAffineField> {
    ground
        .bodies
        .iter()
        .map(|body| MccAffineField {
            body: body.body.clone(),
            reference_position_m: [0.; 3],
            effective_stress_at_reference_compression_positive_kpa: [100., 100., 100., 0., 0., 0.],
            effective_stress_gradient_kpa_per_m: [[0.; 3]; 6],
            pc_over_mean_effective: 1.4,
            void_ratio: VoidRatio::FromMccReference,
        })
        .collect()
}

#[test]
fn affine_seed_matches_original_plastic_history_and_checkpoint_without_mesh_point_values()
-> Result<()> {
    let (mesh, mut ground, case) = fixture(false)?;
    let old = geostatic::states(&mesh, &ground)?;
    let baseline = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(baseline.complete);
    ground.initialization.geostatic.as_mut().unwrap().seed = Seed::MccAffineFields {
        fields: fields(&ground),
    };
    assert_eq!(
        serde_json::to_vec(&old)?,
        serde_json::to_vec(&geostatic::states(&mesh, &ground)?)?
    );
    let run = stages::solve(&mesh, ground.clone(), case.clone())?;
    assert!(run.complete);
    assert_eq!(
        serde_json::to_vec(&baseline.frames)?,
        serde_json::to_vec(&run.frames)?
    );
    assert_eq!(
        serde_json::to_vec(&baseline.last_committed)?,
        serde_json::to_vec(&run.last_committed)?
    );
    let saved = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    assert!(!saved.complete);
    let encoded = serde_json::to_vec(&Checkpoint::seal(&mesh, saved)?)?;
    let cp: Checkpoint = serde_json::from_slice(&encoded)?;
    let resumed = cp.resume(&mesh, ground.clone(), case.clone(), None)?;
    assert_eq!(
        serde_json::to_vec(&run.frames)?,
        serde_json::to_vec(&resumed.frames)?
    );
    assert_eq!(
        serde_json::to_vec(&run.last_committed)?,
        serde_json::to_vec(&resumed.last_committed)?
    );
    let cp: Checkpoint = serde_json::from_slice(&encoded)?;
    if let Seed::MccAffineFields { fields } =
        &mut ground.initialization.geostatic.as_mut().unwrap().seed
    {
        fields[0].effective_stress_gradient_kpa_per_m[2][2] = -0.01;
    }
    assert!(
        cp.resume(&mesh, ground, case, None).is_err(),
        "changed recipe must invalidate saved history"
    );
    Ok(())
}

#[test]
fn affine_global_tensor_uses_all_coordinates_material_ncl_and_rejects_invalid_recipes() -> Result<()>
{
    for quadratic in [false, true] {
        let (mesh, mut ground, _) = fixture(quadratic)?;
        let mut recipe = fields(&ground);
        for f in &mut recipe {
            f.reference_position_m = [1., 2., 3.];
            f.effective_stress_gradient_kpa_per_m = [
                [2., 0., 0.],
                [0., 3., 0.],
                [0., 0., -4.],
                [0.1, 0.2, 0.3],
                [-0.2, 0.3, 0.1],
                [0.3, -0.1, 0.2],
            ];
            f.pc_over_mean_effective = 2.;
        }
        ground.initialization.geostatic.as_mut().unwrap().seed = Seed::MccAffineFields {
            fields: recipe.clone(),
        };
        let state = geostatic::states(&mesh, &ground)?;
        for (s, q) in state.iter().zip(mesh.integration_layout()?) {
            let s = s.as_mcc().unwrap();
            let [x, y, z] = q.position_m;
            let expected = [
                100. + 2. * (x - 1.),
                100. + 3. * (y - 2.),
                100. - 4. * (z - 3.),
                0.1 * (x - 1.) + 0.2 * (y - 2.) + 0.3 * (z - 3.),
                -0.2 * (x - 1.) + 0.3 * (y - 2.) + 0.1 * (z - 3.),
                0.3 * (x - 1.) - 0.1 * (y - 2.) + 0.2 * (z - 3.),
            ];
            for (a, b) in s
                .effective_stress_compression_positive_kpa
                .iter()
                .zip(expected)
            {
                assert!((a - b).abs() < 1e-12);
            }
            let p = expected[..3].iter().sum::<f64>() / 3.;
            assert_eq!(s.preconsolidation_pressure_kpa, 2. * p);
            let m = ground
                .materials
                .iter()
                .find(|m| m.id == mesh.material_ids[q.element_index])
                .unwrap()
                .model
                .as_mcc()
                .unwrap();
            let v = m.specific_volume_ncl_at_reference
                - m.lambda * (2. * p / m.reference_pressure_kpa).ln()
                + m.kappa * 2_f64.ln();
            assert!((s.void_ratio - (v - 1.)).abs() < 1e-14);
            assert_eq!(s.plastic_strain_compression_positive, [0.; 6]);
        }
        for failure in 0..9 {
            let mut bad = recipe.clone();
            match failure {
                0 => {
                    bad.pop();
                }
                1 => bad.push(bad[0].clone()),
                2 => bad[0].body = "foreign".into(),
                3 => bad[0].effective_stress_gradient_kpa_per_m[0][0] = f64::NAN,
                4 => bad[0].pc_over_mean_effective = 0.5,
                5 => bad[0].effective_stress_at_reference_compression_positive_kpa = [-100.; 6],
                6 => bad[0].effective_stress_at_reference_compression_positive_kpa[3] = 1000.,
                7 => bad[0].void_ratio = VoidRatio::Explicit { value: 1. },
                _ => bad[0].reference_position_m[0] = f64::INFINITY,
            }
            ground.initialization.geostatic.as_mut().unwrap().seed =
                Seed::MccAffineFields { fields: bad };
            assert!(
                geostatic::states(&mesh, &ground).is_err(),
                "invalid affine seed {failure} accepted"
            );
        }
        let mut mc: Ground =
            serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
        mc.initialization.geostatic.as_mut().unwrap().seed =
            Seed::MccAffineFields { fields: recipe };
        assert!(geostatic::states(&mesh, &mc).is_err());
    }
    Ok(())
}
