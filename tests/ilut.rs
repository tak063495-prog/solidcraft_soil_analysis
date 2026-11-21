use soil_fem_3d::{
    Result,
    checkpoint::{Checkpoint, hash_json},
    ground::Ground,
    linalg::{CsrMatrix, CsrPattern, GmresOptions, IlutOptions, LinearSolver, gmres_ilu0},
    mesh::Mesh,
    stages::{self, StageCase},
};

#[test]
fn legacy_options_json_has_identical_bytes_when_ilut_is_absent() -> Result<()> {
    let old = r#"{"restart":60,"maximum_iterations":3000,"pivot_tolerance":1e-14}"#;
    let options: GmresOptions = serde_json::from_str(old)?;
    assert!(options.ilut.is_none());
    assert_eq!(serde_json::to_string(&options)?, old);
    assert_eq!(serde_json::to_string(&GmresOptions::default())?, old);
    Ok(())
}

#[test]
fn ilut_rejects_invalid_ranges_types_and_unknown_options() -> Result<()> {
    for drop in [-1., 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            IlutOptions {
                drop_tolerance: drop,
                fill_limit: 20
            }
            .validate()
            .is_err()
        );
    }
    for fill in [0, 4097, usize::MAX] {
        assert!(
            IlutOptions {
                drop_tolerance: 1e-4,
                fill_limit: fill
            }
            .validate()
            .is_err()
        );
    }
    for bad in [
        r#"{"drop_tolerance":0.0001,"fill_limit":1.5}"#,
        r#"{"drop_tolerance":0.0001,"fill_limit":-1}"#,
        r#"{"drop_tolerance":0.0001,"fill_limit":2,"ordering":"colamd"}"#,
    ] {
        assert!(serde_json::from_str::<IlutOptions>(bad).is_err());
    }
    for (drop, fill) in [(0., 1), (1., 4096)] {
        IlutOptions {
            drop_tolerance: drop,
            fill_limit: fill,
        }
        .validate()?;
    }
    Ok(())
}

#[test]
fn ilut_solves_unsymmetric_operator_using_original_residual() -> Result<()> {
    let n = 24;
    let p = CsrPattern::from_connectivity(n, [(0..n).collect()])?;
    let mut a = CsrMatrix::zeros(p);
    for i in 0..n {
        for j in 0..n {
            let value = if i == j {
                30.
            } else {
                ((13 * i + 7 * j) as f64).sin()
            };
            a.add(i, j, value * 10f64.powi((i % 7) as i32 - 3))?;
        }
    }
    let exact: Vec<_> = (0..n).map(|i| (i as f64).cos()).collect();
    let b = a.matvec(&exact)?;
    let options = GmresOptions {
        ilut: Some(IlutOptions {
            drop_tolerance: 0.,
            fill_limit: n,
        }),
        ..Default::default()
    };
    let (x, report) = gmres_ilu0(&a, &b, 1e-10, &options)?;
    assert_eq!(report.method, "csr_gmres_ilut");
    assert!(report.iterations <= 2);
    assert!(a.true_residual(&x, &b)? <= 1e-10);
    Ok(())
}

fn fixture() -> Result<(Mesh, Ground, StageCase)> {
    let mesh = serde_json::from_str(include_str!("../examples/mc-column-mesh.json"))?;
    let ground = serde_json::from_str(include_str!("../examples/ground-geostatic-mc-column.json"))?;
    let mut case: StageCase =
        serde_json::from_str(include_str!("../examples/stages-geostatic-column.json"))?;
    case.settings.linear_solver = LinearSolver::SparseGmresIlu0 {
        options: GmresOptions {
            ilut: Some(IlutOptions {
                drop_tolerance: 1e-4,
                fill_limit: 20,
            }),
            ..Default::default()
        },
    };
    Ok((mesh, ground, case))
}

#[test]
fn checkpoint_binds_ilut_drop_fill_and_preconditioner_choice() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let run = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let cp = Checkpoint::seal(&mesh, run)?;
    cp.verify(&mesh, &ground, &case)?;
    let saved = serde_json::to_value(&cp)?;
    assert_eq!(
        saved["run"]["case"]["settings"]["linear_solver"]["options"]["ilut"]["fill_limit"],
        20
    );
    for change in 0..3 {
        let mut other = case.clone();
        let LinearSolver::SparseGmresIlu0 { options } = &mut other.settings.linear_solver else {
            unreachable!()
        };
        match change {
            0 => options.ilut.as_mut().unwrap().drop_tolerance = 1e-3,
            1 => options.ilut.as_mut().unwrap().fill_limit = 21,
            _ => options.ilut = None,
        }
        assert!(cp.verify(&mesh, &ground, &other).is_err());
    }
    Ok(())
}

#[test]
fn ilut_checkpoint_rejects_rehashed_saved_state_and_saved_settings() -> Result<()> {
    let (mesh, ground, case) = fixture()?;
    let run = stages::solve_until(&mesh, ground.clone(), case.clone(), 1)?;
    let cp = Checkpoint::seal(&mesh, run)?;
    let bytes = serde_json::to_vec(&cp)?;
    for change in 0..2 {
        let mut bad: Checkpoint = serde_json::from_slice(&bytes)?;
        if change == 0 {
            bad.run
                .geostatic
                .as_mut()
                .unwrap()
                .last_committed
                .displacement_m[0][0] += 1.;
        } else {
            let LinearSolver::SparseGmresIlu0 { options } =
                &mut bad.run.case.settings.linear_solver
            else {
                unreachable!()
            };
            options.ilut.as_mut().unwrap().fill_limit = 21;
        }
        bad.run_sha256 = hash_json(&bad.run)?;
        assert!(bad.verify(&mesh, &ground, &case).is_err());
    }
    Ok(())
}
