use soil_fem_3d::{
    Result,
    ground::Ground,
    mc::{Flow, Parameters},
    mesh::Mesh,
    soil,
    stages::{self, McParameterMechanics, StageCase},
};

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mesh: Mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let ground = serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let case = serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    Ok((mesh.elevate_straight()?, ground, case))
}

#[test]
fn unit_epoch_preserves_loaded_plastic_origin_and_independent_controllers() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let origin = stages::solve(&mesh, ground, case)?;
    assert!(origin.complete);
    assert!(
        origin
            .last_committed
            .material_states
            .iter()
            .any(|s| s.plastic_strain().iter().any(|v| v.abs() > 1e-7))
    );
    let saved = serde_json::to_vec(&origin.last_committed)?;
    let mut a = McParameterMechanics::new(&mesh, &origin)?;
    let b = McParameterMechanics::new(&mesh, &origin)?;
    let step = a.try_parameters(a.parameters().to_vec())?;
    assert_eq!(step.iterations, 0);
    assert_eq!(serde_json::to_vec(a.committed())?, saved);
    assert_eq!(serde_json::to_vec(b.committed())?, saved);
    assert_eq!(serde_json::to_vec(&origin.last_committed)?, saved);
    Ok(())
}

fn confined_origin() -> Result<(Mesh, stages::Run)> {
    let (mesh, mut ground, mut case) = fixture()?;
    for m in &mut ground.materials {
        // Positive by input contract; the resulting stress variation is <1e-10 kPa.
        m.unit_weight_kn_m3 = 1e-10;
        m.model = soil::Model::MohrCoulomb(Parameters {
            young_kpa: 30_000.,
            poisson: 0.25,
            cohesion_kpa: 40.,
            friction_deg: 0.,
            dilation_deg: 0.,
            tension_cutoff_kpa: None,
            flow: Flow::NonAssociated,
        });
    }
    case.stages =
        serde_json::from_str(r#"[{"kind":"pressure","name":"load","target_pressure_kpa":100.0}]"#)?;
    let origin = stages::solve(&mesh, ground, case)?;
    assert!(origin.complete, "{:?}", origin.failure_reason);
    Ok((mesh, origin))
}

#[test]
fn post_load_tresca_weakening_matches_confined_compression_closed_form() -> Result<()> {
    let (mesh, origin) = confined_origin()?;
    let mut one = McParameterMechanics::new(&mesh, &origin)?;
    let mut two = McParameterMechanics::new(&mesh, &origin)?;
    let to = |c: f64| {
        one.parameters()
            .iter()
            .cloned()
            .map(|mut p| {
                p.cohesion_kpa = c;
                p
            })
            .collect::<Vec<_>>()
    };
    let intermediate = to(30.);
    let final_parameters = to(20.);
    let step = one.try_parameters(final_parameters.clone())?;
    assert!(step.iterations > 0);
    two.try_parameters(intermediate)?;
    let incremental = two.try_parameters(final_parameters)?;
    // Equilibrium fixes sigma_z=P. Tresca yield fixes sigma_x=sigma_y=P-2c.
    // Associated triaxial plastic flow is isochoric. With epsilon_x=y=0:
    // ep_x=-((1-nu)*sigma_x-nu*P)/E, ep_z=-2*ep_x.
    let p = 100.;
    let c = 20.;
    let e = 30_000.;
    let nu = 0.25;
    let lateral = p - 2. * c;
    let ep_lateral = -((1. - nu) * lateral - nu * p) / e;
    let ep_axial = -2. * ep_lateral;
    let axial = (p - 2. * nu * lateral) / e + ep_axial;
    let original_axial = (p - 2. * nu * (nu * p / (1. - nu))) / e;
    for (i, (coordinate, u)) in mesh
        .nodes
        .iter()
        .zip(&step.frame.state.displacement_m)
        .enumerate()
    {
        assert!(u[0].abs() < 1e-12 && u[1].abs() < 1e-12);
        assert!((u[2] + axial * coordinate[2]).abs() < 1e-9);
        assert!(
            (u[2] - origin.last_committed.displacement_m[i][2]
                + (axial - original_axial) * coordinate[2])
                .abs()
                < 1e-9
        );
        assert!((u[2] - incremental.frame.state.displacement_m[i][2]).abs() < 1e-9);
    }
    for (s, t) in step
        .frame
        .state
        .material_states
        .iter()
        .zip(&incremental.frame.state.material_states)
    {
        let soil::State::MohrCoulomb(s) = s else {
            panic!("expected MC")
        };
        for (actual, expected) in s
            .effective_stress_compression_positive_kpa
            .iter()
            .zip([lateral, lateral, p, 0., 0., 0.])
        {
            assert!((actual - expected).abs() < 1e-5);
        }
        for (actual, expected) in s
            .plastic_strain_compression_positive
            .iter()
            .zip([ep_lateral, ep_lateral, ep_axial, 0., 0., 0.])
        {
            assert!((actual - expected).abs() < 1e-9);
        }
        for (a, b) in s
            .effective_stress_compression_positive_kpa
            .iter()
            .zip(t.stress())
        {
            assert!((a - b).abs() < 1e-5);
        }
    }
    assert_eq!(
        step.frame.state.applied_nodal_load_kn,
        origin.last_committed.applied_nodal_load_kn
    );
    assert_eq!(
        step.frame.state.active_elements,
        origin.last_committed.active_elements
    );
    assert_eq!(
        step.frame.state.stage_index,
        origin.last_committed.stage_index
    );
    Ok(())
}

#[test]
fn invalid_epoch_retains_parameters_and_history_and_rejects_wrong_origin() -> Result<()> {
    let (mesh, mut origin) = confined_origin()?;
    let mut controller = McParameterMechanics::new(&mesh, &origin)?;
    let state = serde_json::to_vec(controller.committed())?;
    let parameters = serde_json::to_vec(controller.parameters())?;
    let mut bad = controller.parameters().to_vec();
    bad[0].young_kpa *= 2.;
    assert!(controller.try_parameters(bad).is_err());
    assert!(controller.try_parameters(vec![]).is_err());
    assert_eq!(serde_json::to_vec(controller.committed())?, state);
    assert_eq!(serde_json::to_vec(controller.parameters())?, parameters);
    drop(controller);
    origin.provenance.mesh_sha256 = "modified".into();
    assert!(McParameterMechanics::new(&mesh, &origin).is_err());
    Ok(())
}
