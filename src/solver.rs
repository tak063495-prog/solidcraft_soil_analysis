use crate::{Error, Result, elastic, mesh::Mesh};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    pub schema: String,
    pub young_kpa: f64,
    pub poisson: f64,
    pub unit_weight_kn_m3: f64,
    pub fixed_node_set: String,
    pub relative_tolerance: f64,
    pub max_iterations: usize,
}

#[derive(Debug, Serialize)]
pub struct Solution {
    pub analysis: &'static str,
    pub iterations: usize,
    pub true_relative_residual: f64,
    pub displacement_m: Vec<[f64; 3]>,
    pub stress_tension_positive_kpa: Vec<[f64; 6]>,
    pub applied_force_kn: [f64; 3],
    pub support_reaction_kn: [f64; 3],
    pub equilibrium_error_kn: [f64; 3],
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn multiply(rows: &[BTreeMap<usize, f64>], x: &[f64]) -> Vec<f64> {
    rows.iter()
        .map(|row| row.iter().map(|(&j, &v)| v * x[j]).sum())
        .collect()
}

pub fn solve(mesh: &Mesh, case: &Case) -> Result<Solution> {
    mesh.validate()?;
    if !mesh.tetrahedra10.is_empty() {
        return Err(Error::Input(
            "bootstrap elastic solver accepts Tet4 only; use solve-mechanics for Tet10".into(),
        ));
    }
    if case.schema != "soil-fem-elastic-case/1"
        || !case.unit_weight_kn_m3.is_finite()
        || case.unit_weight_kn_m3 <= 0.0
        || !case.relative_tolerance.is_finite()
        || !(1e-14..=1e-4).contains(&case.relative_tolerance)
        || case.max_iterations == 0
        || case.max_iterations > 100_000
    {
        return Err(Error::Input("invalid elastic case settings".into()));
    }
    if mesh.material_ids.iter().any(|&id| id != 1) {
        return Err(Error::Input(
            "bootstrap solver supports material 1 only".into(),
        ));
    }
    let fixed = mesh
        .node_sets
        .get(&case.fixed_node_set)
        .ok_or_else(|| Error::Input("fixed node set missing".into()))?;
    if fixed.len() < 3 {
        return Err(Error::Input(
            "at least three non-collinear fixed nodes required".into(),
        ));
    }
    let p0 = nalgebra::Vector3::from_row_slice(&mesh.nodes[fixed[0]]);
    let has_area = fixed.iter().skip(1).any(|&i| {
        fixed.iter().skip(1).any(|&j| {
            let a = nalgebra::Vector3::from_row_slice(&mesh.nodes[i]) - p0;
            let b = nalgebra::Vector3::from_row_slice(&mesh.nodes[j]) - p0;
            a.cross(&b).norm() > 1e-12 * a.norm() * b.norm()
        })
    });
    if !has_area {
        return Err(Error::Input("fixed nodes are collinear".into()));
    }
    let n = mesh.nodes.len() * 3;
    let d = elastic::constitutive(case.young_kpa, case.poisson)?;
    let mut k = vec![BTreeMap::<usize, f64>::new(); n];
    let mut f = vec![0.; n];
    for tet in &mesh.tetrahedra {
        let (v, _, ke) = elastic::element(&mesh.points(tet), &d)?;
        for i in 0..12 {
            let gi = tet[i / 3] * 3 + i % 3;
            if i % 3 == 2 {
                f[gi] -= case.unit_weight_kn_m3 * v / 4.;
            }
            for j in 0..12 {
                *k[gi].entry(tet[j / 3] * 3 + j % 3).or_default() += ke[(i, j)];
            }
        }
    }
    let mut constrained = vec![false; n];
    for &i in fixed {
        constrained[3 * i..3 * i + 3].fill(true);
    }
    let free: Vec<_> = (0..n).filter(|&i| !constrained[i]).collect();
    if free.is_empty() {
        return Err(Error::Input("no free degrees of freedom".into()));
    }
    let mut rhs = f.clone();
    for i in 0..n {
        if constrained[i] {
            rhs[i] = 0.;
        }
    }
    let norm_f = dot(&rhs, &rhs).sqrt();
    if !norm_f.is_finite() || norm_f == 0. {
        return Err(Error::Input("no finite free load".into()));
    }
    let mut diag = vec![1.; n];
    for &i in &free {
        diag[i] = *k[i]
            .get(&i)
            .ok_or_else(|| Error::Solve("missing diagonal".into()))?;
        if !diag[i].is_finite() || diag[i] <= 0. {
            return Err(Error::Solve("invalid diagonal".into()));
        }
    }
    let mut u = vec![0.; n];
    let mut r = rhs.clone();
    let mut z: Vec<_> = (0..n).map(|i| r[i] / diag[i]).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut converged = None;
    for iteration in 1..=case.max_iterations {
        let mut ap = multiply(&k, &p);
        for i in 0..n {
            if constrained[i] {
                ap[i] = 0.;
            }
        }
        let pap = dot(&p, &ap);
        if !pap.is_finite() || pap <= 0. {
            return Err(Error::Solve("CG breakdown; inspect supports/mesh".into()));
        }
        let alpha = rz / pap;
        for &i in &free {
            u[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        if dot(&r, &r).sqrt() / norm_f <= case.relative_tolerance {
            // Acceptance uses the original matrix, not only CG's recursive residual.
            let ku = multiply(&k, &u);
            let true_r = free
                .iter()
                .map(|&i| (rhs[i] - ku[i]).powi(2))
                .sum::<f64>()
                .sqrt()
                / norm_f;
            if true_r.is_finite() && true_r <= case.relative_tolerance {
                converged = Some((iteration, true_r));
                break;
            }
            return Err(Error::Solve(format!("true residual too large: {true_r:e}")));
        }
        for &i in &free {
            z[i] = r[i] / diag[i];
        }
        let next = dot(&r, &z);
        let beta = next / rz;
        for &i in &free {
            p[i] = z[i] + beta * p[i];
        }
        rz = next;
    }
    let (iterations, residual) = converged.ok_or_else(|| Error::Solve("iteration limit".into()))?;
    if !u.iter().all(|x| x.is_finite()) {
        return Err(Error::Solve("non-finite solution".into()));
    }
    let ku = multiply(&k, &u);
    let mut applied = [0.; 3];
    let mut reaction = [0.; 3];
    for i in 0..n {
        applied[i % 3] += f[i];
        if constrained[i] {
            reaction[i % 3] += ku[i] - f[i];
        }
    }
    let stress = mesh
        .tetrahedra
        .iter()
        .map(|tet| {
            let (_, b, _) = elastic::element(&mesh.points(tet), &d)?;
            let ue = std::array::from_fn(|i| u[tet[i / 3] * 3 + i % 3]);
            Ok(elastic::stress(&b, &d, &ue))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Solution {
        analysis: "linear-elastic-tet4-self-weight-bootstrap",
        iterations,
        true_relative_residual: residual,
        displacement_m: u.as_chunks::<3>().0.to_vec(),
        stress_tension_positive_kpa: stress,
        applied_force_kn: applied,
        support_reaction_kn: reaction,
        equilibrium_error_kn: std::array::from_fn(|i| applied[i] + reaction[i]),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn constrained_tetra_has_analytic_settlement_and_reaction() -> Result<()> {
        let mesh = Mesh {
            schema: "soil-fem-mesh/1".into(),
            units: "m-kN-kPa".into(),
            axis: "z_up".into(),
            crs: "LOCAL".into(),
            origin_m: [0.; 3],
            nodes: vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            tetrahedra: vec![[0, 1, 2, 3]],
            tetrahedra10: vec![],
            quadrature_order: None,
            face_sets: BTreeMap::new(),
            material_ids: vec![1],
            node_sets: BTreeMap::from([("base".into(), vec![0, 1, 2])]),
            regions: vec![],
            region_ids: vec![],
            excavation_sets: BTreeMap::new(),
        };
        let case = Case {
            schema: "soil-fem-elastic-case/1".into(),
            young_kpa: 1000.,
            poisson: 0.25,
            unit_weight_kn_m3: 20.,
            fixed_node_set: "base".into(),
            relative_tolerance: 1e-10,
            max_iterations: 100,
        };
        let s = solve(&mesh, &case)?;
        // k_zz=V*(lambda+2mu)=200; f_z=-gamma*V/4=-5/6.
        assert!((s.displacement_m[3][2] + 1. / 240.).abs() < 1e-12);
        assert!((s.support_reaction_kn[2] - 20. / 6.).abs() < 1e-12);
        assert!(s.equilibrium_error_kn.iter().all(|x| x.abs() < 1e-12));
        Ok(())
    }
}
