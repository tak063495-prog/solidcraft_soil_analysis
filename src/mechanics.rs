//! Multi-material Tet4/Tet10 elasticity with general affine constraints and
//! consistent nodal, body, traction and signed pressure loads.
use crate::{
    Error, Result,
    boundary::{Conditions, Transform},
    elastic,
    elements::PointInfo,
    linalg::{self, CsrMatrix, CsrPattern, LinearReport, LinearSolver},
    mesh::Mesh,
};
use nalgebra::{DMatrix, DVector, SVector, Vector3};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    time::Instant,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Material {
    pub id: u32,
    pub young_kpa: f64,
    pub poisson: f64,
    pub body_force_kn_m3: [f64; 3],
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub schema: String,
    pub units: String,
    pub materials: Vec<Material>,
    pub boundary: Conditions,
    pub relative_linear_tolerance: f64,
    pub maximum_free_dofs: usize,
    pub linear_solver: LinearSolver,
}
#[derive(Debug, Serialize)]
pub struct Solution {
    pub schema: &'static str,
    pub solver_binary_sha256: String,
    pub mesh_sha256: String,
    pub case_sha256: String,
    pub case: Case,
    pub displacement_m: Vec<[f64; 3]>,
    pub reaction_kn: Vec<[f64; 3]>,
    pub applied_nodal_load_kn: Vec<[f64; 3]>,
    pub integration_points: Vec<PointInfo>,
    pub stress_tension_positive_kpa: Vec<[f64; 6]>,
    pub strain_engineering: Vec<[f64; 6]>,
    pub true_relative_free_force_residual: f64,
    pub maximum_constraint_error_m: f64,
    pub force_balance_error_kn: [f64; 3],
    pub moment_balance_error_kn_m: [f64; 3],
    pub strain_energy_kn_m: f64,
    pub linear: LinearReport,
}

pub fn solve(mesh: &Mesh, case: Case) -> Result<Solution> {
    mesh.validate()?;
    case.boundary.validate_definition(mesh)?;
    let solver_binary_sha256 = crate::checkpoint::hash_file(&std::env::current_exe()?)?;
    let mesh_sha256 = crate::checkpoint::hash_json(mesh)?;
    let case_sha256 = crate::checkpoint::hash_json(&case)?;
    case.linear_solver.validate()?;
    if case.schema != "soil-fem-mechanics/1"
        || case.units != "m-kN-kPa"
        || !case.relative_linear_tolerance.is_finite()
        || !(1e-12..=1e-5).contains(&case.relative_linear_tolerance)
        || case.maximum_free_dofs == 0
        || case.maximum_free_dofs > case.linear_solver.dimension_limit()
    {
        return Err(Error::Input("invalid general mechanics case".into()));
    }
    let mut materials = BTreeMap::new();
    for m in &case.materials {
        if m.id == 0
            || !m.body_force_kn_m3.iter().all(|v| v.is_finite())
            || materials.insert(m.id, m).is_some()
        {
            return Err(Error::Input("invalid/duplicate mechanics material".into()));
        }
        elastic::constitutive(m.young_kpa, m.poisson)?;
    }
    if mesh.material_ids.iter().copied().collect::<BTreeSet<_>>()
        != materials.keys().copied().collect()
    {
        return Err(Error::Input(
            "mechanics material/mesh assignment mismatch".into(),
        ));
    }
    let transform = Transform::new(
        mesh,
        &case.boundary.constraints,
        vec![true; mesh.tetrahedra.len()],
    )?;
    if transform.count > case.maximum_free_dofs {
        return Err(Error::Input(
            "general mechanics free-DOF limit exceeded".into(),
        ));
    }
    let points: Vec<_> = (0..mesh.tetrahedra.len())
        .map(|i| mesh.integration(i))
        .collect::<Result<_>>()?;
    let mut load = case
        .boundary
        .loads(mesh, &vec![true; mesh.tetrahedra.len()])?;
    let mut lifted_force = vec![0.; load.len()];
    let connectivity: Vec<Vec<usize>> = (0..mesh.tetrahedra.len())
        .map(|i| {
            let mut rows = Vec::new();
            for n in mesh.cell_nodes(i) {
                for a in 0..3 {
                    rows.extend(transform.rows[3 * n + a].iter().map(|&(r, _)| r));
                }
            }
            rows
        })
        .collect();
    let mut sparse = if transform.count > 0 && !case.linear_solver.is_dense() {
        Some(CsrMatrix::zeros(CsrPattern::from_connectivity(
            transform.count,
            connectivity,
        )?))
    } else {
        None
    };
    let mut dense = if transform.count > 0 && case.linear_solver.is_dense() {
        Some(DMatrix::<f64>::zeros(transform.count, transform.count))
    } else {
        None
    };
    for (ei, quadrature) in points.iter().enumerate() {
        let nodes = mesh.cell_nodes(ei);
        let ndof = nodes.len() * 3;
        let dofs: Vec<_> = nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let material = materials[&mesh.material_ids[ei]];
        let d = elastic::constitutive(material.young_kpa, material.poisson)?;
        let mut ke = DMatrix::<f64>::zeros(ndof, ndof);
        for q in quadrature {
            ke += q.weight_m3 * q.b.transpose() * d * &q.b;
            for (i, &node) in nodes.iter().enumerate() {
                for a in 0..3 {
                    load[3 * node + a] += q.weight_m3 * q.shape[i] * material.body_force_kn_m3[a];
                }
            }
        }
        let u0 = DVector::from_iterator(ndof, dofs.iter().map(|&g| transform.offset[g]));
        let f0 = &ke * u0;
        for (i, &gi) in dofs.iter().enumerate() {
            lifted_force[gi] += f0[i];
            for (j, &gj) in dofs.iter().enumerate() {
                for &(ri, ci) in &transform.rows[gi] {
                    for &(rj, cj) in &transform.rows[gj] {
                        if let Some(a) = &mut sparse {
                            a.add(ri, rj, ci * cj * ke[(i, j)])?;
                        }
                        if let Some(a) = &mut dense {
                            a[(ri, rj)] += ci * cj * ke[(i, j)];
                        }
                    }
                }
            }
        }
    }
    let physical_rhs: Vec<_> = load.iter().zip(lifted_force).map(|(f, k)| f - k).collect();
    let rhs = transform.reduce(&physical_rhs);
    let (unknowns, linear) = if let Some(a) = sparse {
        let LinearSolver::SparseGmresIlu0 { options } = &case.linear_solver else {
            return Err(Error::Input("linear backend mismatch".into()));
        };
        linalg::gmres_ilu0(&a, &rhs, case.relative_linear_tolerance, options)?
    } else if let Some(a) = dense {
        let started = Instant::now();
        let b = DVector::from_row_slice(&rhs);
        let x = a
            .clone()
            .lu()
            .solve(&b)
            .ok_or_else(|| Error::Solve("singular constrained mechanics stiffness".into()))?;
        let residual = (&a * &x - b).norm() / DVector::from_row_slice(&rhs).norm().max(1e-30);
        if !residual.is_finite() || residual > case.relative_linear_tolerance {
            return Err(Error::Solve("dense mechanics true-residual failure".into()));
        }
        let r = LinearReport {
            method: "dense_lu".into(),
            dimension: a.nrows(),
            stored_nonzeros: a.len(),
            iterations: 1,
            restarts: 0,
            true_relative_residual: residual,
            matrix_storage_bytes: a.len() * 8,
            workspace_storage_bytes: (a.len() + 2 * a.nrows()) * 8,
            elapsed_seconds: started.elapsed().as_secs_f64(),
        };
        (x.as_slice().to_vec(), r)
    } else {
        (
            vec![],
            LinearReport {
                method: "fully_prescribed".into(),
                dimension: 0,
                stored_nonzeros: 0,
                iterations: 0,
                restarts: 0,
                true_relative_residual: 0.,
                matrix_storage_bytes: 0,
                workspace_storage_bytes: 0,
                elapsed_seconds: 0.,
            },
        )
    };
    let mut u = transform.expand(&unknowns);
    for (v, lift) in u.iter_mut().zip(&transform.offset) {
        *v += lift;
    }
    let mut internal = vec![0.; u.len()];
    let mut stress = Vec::new();
    let mut strain = Vec::new();
    let mut integration_points = Vec::new();
    let mut energy = 0.;
    for (ei, quadrature) in points.iter().enumerate() {
        let nodes = mesh.cell_nodes(ei);
        let dofs: Vec<_> = nodes
            .iter()
            .flat_map(|&n| [3 * n, 3 * n + 1, 3 * n + 2])
            .collect();
        let ue = DVector::from_iterator(dofs.len(), dofs.iter().map(|&g| u[g]));
        let material = materials[&mesh.material_ids[ei]];
        let d = elastic::constitutive(material.young_kpa, material.poisson)?;
        for (qi, q) in quadrature.iter().enumerate() {
            let eps = &q.b * &ue;
            let sig = d * SVector::<f64, 6>::from_column_slice(eps.as_slice());
            let f = q.weight_m3 * q.b.transpose() * sig;
            for (i, &g) in dofs.iter().enumerate() {
                internal[g] += f[i];
            }
            energy += 0.5 * q.weight_m3 * eps.dot(&sig);
            strain.push(std::array::from_fn(|i| eps[i]));
            stress.push(std::array::from_fn(|i| sig[i]));
            integration_points.push(PointInfo {
                element_index: ei,
                local_index: qi,
                natural_coordinates: q.natural,
                position_m: q.position,
                weight_m3: q.weight_m3,
                shape_functions: q.shape.clone(),
            });
        }
    }
    let imbalance: Vec<_> = internal.iter().zip(&load).map(|(i, f)| i - f).collect();
    let reduced = transform.reduce(&imbalance);
    let norm = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let relative = norm(&reduced) / norm(&load).max(norm(&internal)).max(1.);
    let constraint_error = transform.maximum_constraint_error(&u);
    if !relative.is_finite()
        || relative > case.relative_linear_tolerance * 10.
        || constraint_error > 1e-10
    {
        return Err(Error::Solve(
            "mechanics force/constraint verification failed".into(),
        ));
    }
    let reactions = transform.reactions(&imbalance);
    let mut force = [0.; 3];
    let mut moment = Vector3::zeros();
    for (node, p) in mesh.nodes.iter().enumerate() {
        let total = Vector3::from_fn(|a, _| load[3 * node + a] + reactions[3 * node + a]);
        for a in 0..3 {
            force[a] += total[a];
        }
        moment += Vector3::from_row_slice(p).cross(&total);
    }
    let packed = |v: &[f64]| v.as_chunks::<3>().0.to_vec();
    Ok(Solution {
        solver_binary_sha256,
        mesh_sha256,
        case_sha256,
        schema: "soil-fem-mechanics-result/1",
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
    })
}
