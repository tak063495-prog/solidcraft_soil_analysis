//! Diagnostic-only bridge to the existing predictor, in an isolated source copy.
use crate::{Result, joints, mesh::Mesh, rigidity};
use nalgebra::{DMatrix, DVector};
use serde_json::{Value, json};

pub fn predictor(v: Value) -> Result<Value> {
    let mesh: Mesh = serde_json::from_value(v["mesh"].clone())?;
    assert!(mesh.nodes.len() * 3 <= 256 && mesh.tetrahedra.len() <= 48);
    let model = joints::Model::new(&mesh, serde_json::from_value(v["configuration"].clone())?)?;
    let active = vec![true; mesh.tetrahedra.len()];
    let rows: Vec<Vec<(usize, f64)>> = serde_json::from_value(v["rows"].clone())?;
    let a: Vec<Vec<f64>> = serde_json::from_value(v["matrix"].clone())?;
    let n = a.len();
    let matrix = DMatrix::from_fn(n, n, |i, j| a[i][j]);
    let rhs = DVector::from_vec(serde_json::from_value(v["rhs"].clone())?);
    let geometry = rigidity::Geometry::new(&mesh, &active)?;
    let mut old = model.initial_state()?;
    let mut current = old.clone();
    for (key, state) in [("old", &mut old), ("current", &mut current)] {
        for (i, life) in state.lifecycle.iter_mut().enumerate() {
            life.active = v["born"].as_bool().unwrap();
            life.tied = !life.active;
            life.generation = v["generations"][i].as_u64().unwrap();
        }
        for (i, point) in state.points.iter_mut().enumerate() {
            point.jump_m = serde_json::from_value(v[key]["jump"][i].clone())?;
            point.traction_kpa = serde_json::from_value(v[key]["traction"][i].clone())?;
            point.plastic_slip_m = serde_json::from_value(v["old"]["plastic"][i].clone())?;
            point.accumulated_plastic_slip_m = v["old"]["accumulated"][i].as_f64().unwrap();
            point.plastic_dissipation_kpa_m = v["old"]["pd"][i].as_f64().unwrap();
            point.numerical_dissipation_kpa_m = v["old"]["nd"][i].as_f64().unwrap();
            point.endpoint_work_kpa_m = v["old"]["work"][i].as_f64().unwrap();
            point.status = match v[key]["status"][i].as_i64().unwrap() {
                0 => crate::joint_law::Status::Open,
                2 => crate::joint_law::Status::Slip,
                _ => crate::joint_law::Status::Stick,
            };
            point.birth = Some(crate::joint_law::BirthReference {
                jump_m: serde_json::from_value(v["old"]["reference"][i].clone())?,
                initial_traction_kpa: serde_json::from_value(v["old"]["initial"][i].clone())?,
                elastic_offset_m: serde_json::from_value(v[key]["rest"][i].clone())?,
                initial_energy_kpa_m: v["old"]["initial_energy"][i].as_f64().unwrap(),
            });
        }
    }
    let trust = !geometry.constrained(&rows)?.full();
    let mut space = geometry.constrained(&rows)?;
    model.add_rigid_support(
        &current,
        &active,
        v["pressure_tolerance"].as_f64().unwrap(),
        &mut space,
    )?;
    if let Some(k) = v["support_k"].as_array() {
        let scale = k
            .iter()
            .flat_map(|r| r.as_array().unwrap())
            .map(|x| x.as_f64().unwrap().abs())
            .fold(0., f64::max);
        if scale > 0. {
            for (i, row) in k.iter().enumerate() {
                let terms: Vec<_> = row
                    .as_array()
                    .unwrap()
                    .iter()
                    .enumerate()
                    .map(|(j, x)| (j, x.as_f64().unwrap() / scale))
                    .collect();
                let _ = i;
                space.physical_row(&terms)?;
            }
        }
    }
    let gauges = if space.full() {
        None
    } else {
        let modes = geometry.reduced_modes(&space.null_modes(), &rows, n)?;
        Some(rigidity::gauge_coordinates(
            &modes,
            |x| Ok(&matrix * x),
            (0..n)
                .map(|i| matrix.row(i).iter().map(|x| x.abs()).sum::<f64>())
                .fold(0., f64::max),
        )?)
    };
    let cold = if v["iteration"] == 0 && gauges.is_some() {
        model.cold_activation_tangent(&old, &active)?
    } else {
        None
    };
    let mut used = false;
    let (x, coordinates) = if let Some(elastic) = cold {
        used = true;
        let u: Vec<[f64; 3]> = serde_json::from_value(v["displacement"].clone())?;
        let counts: Vec<usize> = serde_json::from_value(v["counts"].clone())?;
        let actual = model.evaluate(&old, &u, &active, Some(&counts), true)?;
        let mut temporary = matrix.clone();
        for (sign, terms) in [(-1., actual.tangent.as_slice()), (1., elastic.as_slice())] {
            for &(g, h, value) in terms {
                for &(i, ci) in &rows[g] {
                    for &(j, cj) in &rows[h] {
                        temporary[(i, j)] += sign * ci * cj * value;
                    }
                }
            }
        }
        (temporary.lu().solve(&rhs), vec![])
    } else if gauges.as_ref().is_none_or(|g| g.coordinates.is_empty()) {
        (matrix.clone().lu().solve(&rhs), vec![])
    } else {
        used = true;
        let g = gauges.as_ref().unwrap();
        let mut b = if v["projected_load"].as_bool().unwrap() {
            &rhs - &g.modes * (g.modes.transpose() * &rhs)
        } else {
            rhs.clone()
        };
        let mut temporary = matrix.clone();
        let scale = (0..n)
            .map(|i| matrix.row(i).iter().map(|x| x.abs()).sum::<f64>())
            .fold(0., f64::max)
            .max(1.);
        for &i in &g.coordinates {
            temporary.row_mut(i).fill(0.);
            temporary.column_mut(i).fill(0.);
            temporary[(i, i)] = scale;
            b[i] = 0.;
        }
        let mut x = temporary
            .lu()
            .solve(&b)
            .ok_or_else(|| crate::Error::Solve("diagnostic predictor LU singular".into()))?;
        x -= &g.modes * (g.modes.transpose() * &x);
        let eta = model.neutral_elastic_predictor(
            (&old, &current),
            &active,
            &rows,
            &g.modes,
            &x,
            g.modes.transpose() * &rhs,
        )?;
        x += &g.modes * eta;
        (Some(x), g.coordinates.clone())
    };
    Ok(match x {
        Some(x) => json!({"correction":x.as_slice(),"used":used,"coordinates":coordinates,
        "trust":trust,"original_linear_relative":(&matrix*&x-&rhs).norm()/rhs.norm().max(1e-30)}),
        None => json!({"error":"diagnostic predictor singular","used":used}),
    })
}
