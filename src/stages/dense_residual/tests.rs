use super::*;
use crate::stages::{Settings, StageCase, Tangent};
use serde::Deserialize;
use serde_json::json;
use std::{fs, path::Path};

const METHODS: [DenseResidualMethod; 4] = [
    DenseResidualMethod::RepeatedRefinement,
    DenseResidualMethod::ExtendedResidual,
    DenseResidualMethod::Ruiz,
    DenseResidualMethod::GmresPolish,
];

fn case() -> Result<StageCase> {
    Ok(serde_json::from_str(include_str!(
        "../../../examples/stages-free-contact.json"
    ))?)
}

#[test]
fn rejects_conflicting_and_unknown_settings() -> Result<()> {
    for method in METHODS
        .into_iter()
        .chain([DenseResidualMethod::CoordinatePolish])
    {
        let mut c = case()?;
        c.settings.dense_residual_method = Some(method);
        assert!(
            matches!(c.validate(), Err(Error::Input(reason)) if reason.contains("requires dense_lu"))
        );
        c.settings.linear_solver = Default::default();
        c.settings.maximum_free_dofs = 750;
        c.validate()?;
        c.settings.dense_lu_refinement = true;
        assert!(
            matches!(c.validate(), Err(Error::Input(reason)) if reason.contains("select only one"))
        );
    }
    let mut v = serde_json::to_value(case()?.settings)?;
    assert!(v.get("dense_residual_method").is_none());
    v["dense_residual_method"] = json!("other");
    assert!(serde_json::from_value::<Settings>(v).is_err());
    Ok(())
}

#[test]
fn recovers_a_low_component_after_two_levels_of_cancellation() {
    let p = 2_f64.powi(120);
    let q = 2_f64.powi(60);
    let a = DMatrix::from_row_slice(1, 5, &[p, q, 1., -q, -p]);
    let x = DVector::from_element(5, 1.);
    let b = DVector::zeros(1);
    assert_eq!(extended_residual(&a, &x, &b)[0], -1.);
    assert_ne!(dense_refinement::compensated_residual(&a, &x, &b)[0], -1.);
}

#[test]
fn all_methods_keep_exact_solutions_and_reject_singular_systems() -> Result<()> {
    let a = DMatrix::from_row_slice(2, 2, &[2., 3., 0., 4.]);
    let b = DVector::from_vec(vec![8., 8.]);
    for method in METHODS
        .into_iter()
        .chain([DenseResidualMethod::CoordinatePolish])
    {
        let (x, report) = solve(&a, &b, 1e-10, method)?;
        assert_eq!(x.as_slice(), &[1., 2.]);
        assert_eq!(report.iterations, 1);
        assert_eq!(report.true_relative_residual, 0.);
        assert!(solve(&DMatrix::zeros(2, 2), &b, 1e-10, method).is_err());
    }
    Ok(())
}

#[derive(Deserialize)]
struct Saved {
    matrix: Vec<Vec<f64>>,
    rhs: Vec<f64>,
}

#[test]
fn replay_all_eight_systems_through_stage_solver() -> Result<()> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut seen = 0;
    let mut settings = case()?.settings;
    settings.linear_solver = Default::default();
    for file in [
        "tests/fixtures/linear-refinement/nc96.jsonl",
        "interim/2026-11-22-impl-linear-solver/refined-terminal.jsonl",
    ] {
        let raw = fs::read(root.join(file))?;
        for (idx, data) in raw
            .split(|&v| v == b'\n')
            .filter(|line| !line.is_empty())
            .enumerate()
        {
            let saved: Saved = serde_json::from_slice(data)?;
            let b = DVector::from_vec(saved.rhs);
            let matrix = Tangent::Dense(DMatrix::from_fn(b.len(), b.len(), |i, j| {
                saved.matrix[i][j]
            }));
            for (name, method) in [("baseline", None), ("refined", None)]
                .into_iter()
                .chain(METHODS.map(|m| (m.name(), Some(m))))
            {
                settings.dense_residual_method = method;
                settings.dense_lu_refinement = name == "refined";
                let start = Instant::now();
                let result = matrix.solve(&b, &settings);
                let seconds = start.elapsed().as_secs_f64();
                let observation = match result {
                    Ok((x, report)) => {
                        let actual = (matrix.action(&x)? - &b).norm() / b.norm().max(1e-30);
                        assert!(x.iter().all(|v| v.is_finite()));
                        assert_eq!(actual, report.true_relative_residual);
                        assert!(actual <= 1e-10);
                        assert_ne!(name, "baseline");
                        json!({"passed":true,"residual":actual,"iterations":report.iterations})
                    }
                    Err(Error::Solve(reason)) => {
                        let residual: f64 = reason
                            .strip_prefix("true linear residual ")
                            .unwrap()
                            .split_whitespace()
                            .next()
                            .unwrap()
                            .parse()
                            .unwrap();
                        assert!(residual > 1e-10);
                        if file.ends_with("nc96.jsonl") {
                            assert!(!matches!(name, "refined" | "repeated_refinement"));
                        }
                        json!({"passed":false,"residual":residual,"reason":reason})
                    }
                    Err(e) => return Err(e),
                };
                println!(
                    "NC01_REPLAY {}",
                    json!({"file":file,"index":idx,"method":name,"seconds":seconds,"result":observation})
                );
            }
            seen += 1;
        }
    }
    assert_eq!(seen, 8);
    Ok(())
}

#[test]
fn coordinate_polish_replays_seven_rejections_and_refined_terminal() -> Result<()> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut settings = case()?.settings;
    settings.linear_solver = Default::default();
    settings.dense_residual_method = Some(DenseResidualMethod::CoordinatePolish);
    assert_eq!(settings.relative_linear_tolerance, 1e-10);
    for (file, count) in [
        ("tests/fixtures/linear-refinement/nc96.jsonl", 7),
        (
            "interim/2026-11-22-impl-linear-solver/refined-terminal.jsonl",
            1,
        ),
    ] {
        let raw = fs::read(root.join(file))?;
        let mut seen = 0;
        for entry in serde_json::Deserializer::from_slice(&raw).into_iter::<Saved>() {
            let saved = entry?;
            let b = DVector::from_vec(saved.rhs);
            let matrix = Tangent::Dense(DMatrix::from_fn(b.len(), b.len(), |i, j| {
                saved.matrix[i][j]
            }));
            let (x, report) = matrix.solve(&b, &settings)?;
            let absolute = (matrix.action(&x)? - &b).norm();
            let relative = absolute / b.norm().max(1e-30);
            assert!(x.iter().all(|v| v.is_finite()));
            assert_eq!(report.true_absolute_residual, Some(absolute));
            assert_eq!(report.true_relative_residual, relative);
            assert!(relative <= settings.relative_linear_tolerance);
            assert!(report.iterations >= 2);
            if count == 1 {
                assert!(report.iterations > 2);
            }
            println!(
                "NC01_COORDINATE {}",
                json!({"file":file,"index":seen,"dimension":b.len(),"absolute":absolute,"relative":relative,"rhs_norm":b.norm(),"iterations":report.iterations})
            );
            seen += 1;
        }
        assert_eq!(seen, count);
    }
    Ok(())
}

#[test]
fn coordinate_sweep_projects_the_residual_along_each_column() {
    let a = DMatrix::from_row_slice(2, 2, &[0., 4., -2., 0.]);
    let b = DVector::from_vec(vec![12., 2.]);
    let system = System {
        a: &a,
        b: &b,
        scale: b.norm(),
        tolerance: 1e-10,
    };
    let mut x = DVector::zeros(2);
    assert_eq!(coordinate::polish(&system, &mut x), (1, 0., 0.));
    assert_eq!(x.as_slice(), &[-1., 3.]);

    let a = DMatrix::from_row_slice(2, 2, &[1., 0., 0., 0.]);
    let system = System { a: &a, ..system };
    let (_, absolute, relative) = coordinate::polish(&system, &mut x);
    assert_eq!(absolute, 2.);
    assert!(relative > system.tolerance);
}
