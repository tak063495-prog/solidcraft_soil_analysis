//! Small-strain Tet10 displacement / corner-P1 mean-pressure elasticity.
//!
//! sigma = 2 G dev(epsilon) - p I, div(u) + p/K = 0. Pressure is
//! compression-positive and split at material interfaces. No penalty shift,
//! reduced integration, stress projection or automatic solver fallback is used.
use crate::{
    Error, Result,
    boundary::Transform,
    elastic,
    elements::PointInfo,
    linalg::{self, CsrMatrix, CsrPattern, LinearReport, LinearSolver},
    mechanics::{Case, Solution as MechanicsSolution},
    mesh::Mesh,
};
use nalgebra::{DMatrix, DVector, SMatrix, SVector, Vector3};
use serde::Serialize;
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

#[derive(Debug, Serialize)]
pub struct PressureDof {
    pub node_index: usize,
    pub material_id: u32,
    pub compression_positive_kpa: f64,
}

#[derive(Debug, Serialize)]
pub struct Solution {
    #[serde(flatten)]
    pub mechanics: MechanicsSolution,
    pub formulation: &'static str,
    pub pressure_dofs: Vec<PressureDof>,
    pub pressure_compression_positive_kpa: Vec<f64>,
    pub pressure_equation_volume_residual_l2_m3: f64,
    pub pressure_equation_relative_residual: f64,
}

fn moduli(young: f64, poisson: f64) -> (f64, f64, SMatrix<f64, 6, 6>) {
    let shear = young / (2. * (1. + poisson));
    let bulk = young / (3. * (1. - 2. * poisson));
    let mut dev = SMatrix::<f64, 6, 6>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            dev[(i, j)] = -2. * shear / 3.;
        }
        dev[(i, i)] += 2. * shear;
        dev[(i + 3, i + 3)] = shear;
    }
    (shear, bulk, dev)
}

pub fn solve(mesh: &Mesh, case: Case) -> Result<Solution> {
    mesh.validate()?;
    case.boundary.validate_definition(mesh)?;
    case.linear_solver.validate()?;
    if mesh.tetrahedra10.len() != mesh.tetrahedra.len() {
        return Err(Error::Input(
            "mixed elasticity requires Tet10 displacement cells".into(),
        ));
    }
    if mesh.quadrature_order.unwrap_or(3) < 3 {
        return Err(Error::Input(
            "mixed elasticity requires quadrature order at least 3".into(),
        ));
    }
    if case.schema != "soil-fem-mechanics/1"
        || case.units != "m-kN-kPa"
        || !case.relative_linear_tolerance.is_finite()
        || !(1e-12..=1e-5).contains(&case.relative_linear_tolerance)
        || case.maximum_free_dofs == 0
        || case.maximum_free_dofs > case.linear_solver.dimension_limit()
    {
        return Err(Error::Input("invalid mixed mechanics case".into()));
    }
    let mut materials = BTreeMap::new();
    for material in &case.materials {
        if material.id == 0
            || !material.body_force_kn_m3.iter().all(|v| v.is_finite())
            || materials.insert(material.id, material).is_some()
        {
            return Err(Error::Input(
                "invalid/duplicate mixed mechanics material".into(),
            ));
        }
        elastic::constitutive(material.young_kpa, material.poisson)?;
    }
    if mesh.material_ids.iter().copied().collect::<BTreeSet<_>>()
        != materials.keys().copied().collect()
    {
        return Err(Error::Input(
            "mixed mechanics material/mesh assignment mismatch".into(),
        ));
    }
    let transform = Transform::new(
        mesh,
        &case.boundary.constraints,
        vec![true; mesh.tetrahedra.len()],
    )?;
    let keys: BTreeSet<_> = mesh
        .tetrahedra
        .iter()
        .zip(&mesh.material_ids)
        .flat_map(|(cell, &material)| cell.iter().map(move |&node| (node, material)))
        .collect();
    let pressure_index: BTreeMap<_, _> = keys
        .into_iter()
        .enumerate()
        .map(|(i, key)| (key, i))
        .collect();
    let np = pressure_index.len();
    let total = transform.count + np;
    if total > case.maximum_free_dofs {
        return Err(Error::Input(
            "mixed displacement/pressure DOF limit exceeded".into(),
        ));
    }
    let connectivity: Vec<_> = mesh
        .tetrahedra10
        .iter()
        .enumerate()
        .map(|(element, nodes)| {
            let mut rows = Vec::new();
            for &node in nodes {
                for a in 0..3 {
                    rows.extend(transform.rows[3 * node + a].iter().map(|&(r, _)| r));
                }
            }
            rows.extend(mesh.tetrahedra[element].iter().map(|&node| {
                transform.count + pressure_index[&(node, mesh.material_ids[element])]
            }));
            rows
        })
        .collect();
    let mut sparse = if case.linear_solver.is_dense() {
        None
    } else {
        Some(CsrMatrix::zeros(CsrPattern::from_connectivity(
            total,
            connectivity,
        )?))
    };
    let mut dense = if case.linear_solver.is_dense() {
        Some(DMatrix::<f64>::zeros(total, total))
    } else {
        None
    };
    let quadrature = (0..mesh.tetrahedra.len())
        .map(|i| mesh.integration(i))
        .collect::<Result<Vec<_>>>()?;
    let mut load = case
        .boundary
        .loads(mesh, &vec![true; mesh.tetrahedra.len()])?;
    let mut lifting = vec![0.; transform.offset.len()];
    let mut pressure_rhs = vec![0.; np];
    for (element, points) in quadrature.iter().enumerate() {
        let nodes = &mesh.tetrahedra10[element];
        let dofs: Vec<_> = nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let material = materials[&mesh.material_ids[element]];
        let (_, bulk, dev) = moduli(material.young_kpa, material.poisson);
        let pi = mesh.tetrahedra[element].map(|node| pressure_index[&(node, material.id)]);
        let mut local = DMatrix::<f64>::zeros(34, 34);
        for point in points {
            let l = [
                1. - point.natural.iter().sum::<f64>(),
                point.natural[0],
                point.natural[1],
                point.natural[2],
            ];
            let stiffness = point.weight_m3 * point.b.transpose() * dev * &point.b;
            for i in 0..30 {
                for j in 0..30 {
                    local[(i, j)] += stiffness[(i, j)];
                }
                let div = point.b[(0, i)] + point.b[(1, i)] + point.b[(2, i)];
                for a in 0..4 {
                    let value = -point.weight_m3 * l[a] * div;
                    local[(i, 30 + a)] += value;
                    local[(30 + a, i)] += value;
                }
            }
            for a in 0..4 {
                for b in 0..4 {
                    local[(30 + a, 30 + b)] -= point.weight_m3 * l[a] * l[b] / bulk;
                }
            }
            for (i, &node) in nodes.iter().enumerate() {
                for a in 0..3 {
                    load[3 * node + a] +=
                        point.weight_m3 * point.shape[i] * material.body_force_kn_m3[a];
                }
            }
        }
        let mut mapping: Vec<_> = dofs.iter().map(|&g| transform.rows[g].clone()).collect();
        mapping.extend(pi.iter().map(|&p| vec![(transform.count + p, 1.)]));
        for i in 0..34 {
            for j in 0..34 {
                if j < 30 {
                    let value = local[(i, j)] * transform.offset[dofs[j]];
                    if i < 30 {
                        lifting[dofs[i]] += value;
                    } else {
                        pressure_rhs[pi[i - 30]] -= value;
                    }
                }
                for &(row, cr) in &mapping[i] {
                    for &(column, cc) in &mapping[j] {
                        let value = cr * cc * local[(i, j)];
                        if let Some(a) = &mut sparse {
                            a.add(row, column, value)?;
                        }
                        if let Some(a) = &mut dense {
                            a[(row, column)] += value;
                        }
                    }
                }
            }
        }
    }
    let physical_rhs: Vec<_> = load.iter().zip(&lifting).map(|(f, k)| f - k).collect();
    let mut rhs = transform.reduce(&physical_rhs);
    rhs.extend(pressure_rhs);
    let (unknowns, linear) = if let Some(a) = sparse {
        let LinearSolver::SparseGmresIlu0 { options } = &case.linear_solver else {
            return Err(Error::Input("mixed sparse backend mismatch".into()));
        };
        linalg::gmres_ilu0(&a, &rhs, case.relative_linear_tolerance, options)?
    } else if let Some(a) = dense {
        let start = Instant::now();
        let b = DVector::from_row_slice(&rhs);
        let x = a
            .clone()
            .lu()
            .solve(&b)
            .ok_or_else(|| Error::Solve("singular mixed elasticity system".into()))?;
        let residual = (&a * &x - &b).norm() / b.norm().max(1e-30);
        if !residual.is_finite() || residual > case.relative_linear_tolerance {
            return Err(Error::Solve("mixed dense true-residual failure".into()));
        }
        let report = LinearReport {
            method: "dense_lu_mixed".into(),
            dimension: total,
            stored_nonzeros: a.len(),
            iterations: 1,
            restarts: 0,
            true_relative_residual: residual,
            matrix_storage_bytes: a.len() * 8,
            workspace_storage_bytes: (a.len() + 2 * total) * 8,
            elapsed_seconds: start.elapsed().as_secs_f64(),
        };
        (x.as_slice().to_vec(), report)
    } else {
        return Err(Error::Solve("missing mixed matrix".into()));
    };
    let p = &unknowns[transform.count..];
    let mut u = transform.expand(&unknowns[..transform.count]);
    for (value, lift) in u.iter_mut().zip(&transform.offset) {
        *value += lift;
    }
    let mut internal = vec![0.; u.len()];
    let mut pressure_equation = vec![0.; np];
    let mut pressure_scale = vec![0.; np];
    let mut integration_points = Vec::new();
    let mut strain = Vec::new();
    let mut stress = Vec::new();
    let mut gp_pressure = Vec::new();
    let mut energy = 0.;
    for (element, points) in quadrature.iter().enumerate() {
        let nodes = &mesh.tetrahedra10[element];
        let dofs: Vec<_> = nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let ue = DVector::from_iterator(30, dofs.iter().map(|&g| u[g]));
        let material = materials[&mesh.material_ids[element]];
        let (_, bulk, dev) = moduli(material.young_kpa, material.poisson);
        let pi = mesh.tetrahedra[element].map(|node| pressure_index[&(node, material.id)]);
        for (qi, point) in points.iter().enumerate() {
            let l = [
                1. - point.natural.iter().sum::<f64>(),
                point.natural[0],
                point.natural[1],
                point.natural[2],
            ];
            let pressure: f64 = (0..4).map(|a| l[a] * p[pi[a]]).sum();
            let eps = &point.b * &ue;
            let mut sig = dev * SVector::<f64, 6>::from_column_slice(eps.as_slice());
            for a in 0..3 {
                sig[a] -= pressure;
            }
            let force = point.weight_m3 * point.b.transpose() * sig;
            for (i, &g) in dofs.iter().enumerate() {
                internal[g] += force[i];
            }
            let div: f64 = eps.iter().take(3).sum();
            for a in 0..4 {
                pressure_equation[pi[a]] += point.weight_m3 * l[a] * (div + pressure / bulk);
                pressure_scale[pi[a]] +=
                    point.weight_m3 * l[a] * (eps.norm() + (pressure / bulk).abs());
            }
            energy += 0.5 * point.weight_m3 * eps.dot(&sig);
            strain.push(std::array::from_fn(|a| eps[a]));
            stress.push(std::array::from_fn(|a| sig[a]));
            gp_pressure.push(pressure);
            integration_points.push(PointInfo {
                element_index: element,
                local_index: qi,
                natural_coordinates: point.natural,
                position_m: point.position,
                weight_m3: point.weight_m3,
                shape_functions: point.shape.clone(),
            });
        }
    }
    let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
    let imbalance: Vec<_> = internal.iter().zip(&load).map(|(i, f)| i - f).collect();
    let relative = norm(&transform.reduce(&imbalance)) / norm(&load).max(norm(&internal)).max(1.);
    let constraint_error = transform.maximum_constraint_error(&u);
    let pressure_absolute = norm(&pressure_equation);
    let volume: f64 = quadrature.iter().flatten().map(|q| q.weight_m3).sum();
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for xyz in &mesh.nodes {
        for a in 0..3 {
            low[a] = low[a].min(xyz[a]);
            high[a] = high[a].max(xyz[a]);
        }
    }
    let length = (0..3)
        .map(|a| (high[a] - low[a]).powi(2))
        .sum::<f64>()
        .sqrt();
    // The displacement scale also covers prescribed rigid motions, whose
    // physical strain is zero up to floating-point cancellation.
    let displacement_scale = volume * norm(&u) / (u.len() as f64).sqrt() / length;
    let pressure_relative =
        pressure_absolute / norm(&pressure_scale).max(displacement_scale).max(1e-30);
    if !relative.is_finite()
        || relative > 10. * case.relative_linear_tolerance
        || constraint_error > 1e-10
        || !pressure_relative.is_finite()
        || pressure_relative > 10. * case.relative_linear_tolerance
    {
        return Err(Error::Solve(
            "mixed force/constraint/pressure verification failed".into(),
        ));
    }
    let reactions = transform.reactions(&imbalance);
    let mut force = [0.; 3];
    let mut moment = Vector3::zeros();
    for (node, xyz) in mesh.nodes.iter().enumerate() {
        let total = Vector3::from_fn(|a, _| load[3 * node + a] + reactions[3 * node + a]);
        for a in 0..3 {
            force[a] += total[a];
        }
        moment += Vector3::from_row_slice(xyz).cross(&total);
    }
    let packed = |v: &[f64]| v.as_chunks::<3>().0.to_vec();
    let pressure_dofs = pressure_index
        .into_iter()
        .map(|((node, material), i)| PressureDof {
            node_index: node,
            material_id: material,
            compression_positive_kpa: p[i],
        })
        .collect();
    Ok(Solution {
        formulation: "Tet10 displacement / material-split corner P1 mean pressure; finite bulk modulus",
        pressure_dofs,
        pressure_compression_positive_kpa: gp_pressure,
        pressure_equation_volume_residual_l2_m3: pressure_absolute,
        pressure_equation_relative_residual: pressure_relative,
        mechanics: MechanicsSolution {
            schema: "soil-fem-mixed-mechanics-result/1",
            solver_binary_sha256: crate::checkpoint::hash_file(&std::env::current_exe()?)?,
            mesh_sha256: crate::checkpoint::hash_json(mesh)?,
            case_sha256: crate::checkpoint::hash_json(&case)?,
            case,
            displacement_m: packed(&u),
            reaction_kn: packed(&reactions),
            applied_nodal_load_kn: packed(&load),
            integration_points,
            stress_tension_positive_kpa: stress,
            strain_engineering: strain,
            true_relative_free_force_residual: relative,
            maximum_constraint_error_m: constraint_error,
            force_balance_error_kn: force,
            moment_balance_error_kn_m: moment.into(),
            strain_energy_kn_m: energy,
            linear,
        },
    })
}
