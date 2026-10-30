use super::*;
use mixed_element::{Element, Volume};
use soil_fem_3d::{elements, soil};
fn element(volume: Volume, curved: bool) -> Result<Element> {
    let mut nodes = vec![
        [0., 0., 0.],
        [1.3, 0.2, 0.1],
        [0.1, 1.1, 0.],
        [0., 0.2, 1.2],
    ];
    for [i, j] in elements::EDGES {
        nodes.push(std::array::from_fn(|a| (nodes[i][a] + nodes[j][a]) / 2.));
    }
    if curved {
        nodes[4][2] += 0.04;
    }
    Element::new(nodes, volume)
}
fn fixture(e: &Element, p: &Parameters) -> Result<(Vec<soil::State>, Vec<f64>)> {
    let old = initial();
    let de = increment(p, &old, [300., 100., 50., 0., 0., 0.]);
    let (dev, v) = split(de);
    let mean = mc_split::integrate(p, &old, dev, v, 1)?.mean_pressure_compression_positive_kpa;
    let tensor = [
        [de[0], de[3] / 2., de[5] / 2.],
        [de[3] / 2., de[1], de[4] / 2.],
        [de[5] / 2., de[4] / 2., de[2]],
    ];
    let mut x = vec![0.; e.count()];
    for value in &mut x[..e.nt] {
        *value = v;
    }
    for (i, node) in e.nodes.iter().enumerate() {
        for a in 0..3 {
            x[e.nt + 3 * i + a] = -(0..3).map(|b| tensor[a][b] * node[b]).sum::<f64>();
        }
    }
    for value in &mut x[e.nt + 30..] {
        *value = mean;
    }
    Ok((vec![soil::State::MohrCoulomb(old); e.points.len()], x))
}
fn check_jacobian(volume: Volume) -> Result<()> {
    let e = element(volume, false)?;
    let p = parameters(5.);
    let model = soil::Model::MohrCoulomb(p.clone());
    let (old, x) = fixture(&e, &p)?;
    let before: Vec<_> = old
        .iter()
        .map(|v| {
            let soil::State::MohrCoulomb(s) = v else {
                panic!()
            };
            bits(s)
        })
        .collect();
    let evaluated = e.evaluate(&model, &old, &x, true)?;
    let matrix = evaluated.tangent.unwrap();
    assert!(
        evaluated
            .states
            .iter()
            .all(|s| matches!(s,soil::State::MohrCoulomb(v) if v.last_return==Return::ShearFace))
    );
    let mut numeric = vec![vec![0.; e.count()]; e.count()];
    for j in 0..e.count() {
        let h = 2e-8 * x[j].abs().max(1.);
        let mut plus = x.clone();
        let mut minus = x.clone();
        plus[j] += h;
        minus[j] -= h;
        let a = e.evaluate(&model, &old, &plus, false)?;
        let b = e.evaluate(&model, &old, &minus, false)?;
        for i in 0..e.count() {
            numeric[i][j] = (a.residual[i] - b.residual[i]) / (2. * h);
        }
    }
    let ranges = [0..e.nt, e.nt..e.nt + 30, e.nt + 30..e.count()];
    for rows in &ranges {
        for cols in &ranges {
            let mut diff = 0.;
            let mut norm = 0.;
            for i in rows.clone() {
                for j in cols.clone() {
                    diff += (matrix[i][j] - numeric[i][j]).powi(2);
                    norm += numeric[i][j].powi(2);
                }
            }
            assert!(
                diff.sqrt() / norm.sqrt().max(1e-10) < 3e-5,
                "Jacobian block {rows:?}/{cols:?}: {}",
                diff.sqrt() / norm.sqrt().max(1e-10)
            );
        }
    }
    let asym = (0..e.count())
        .flat_map(|i| (0..e.count()).map(move |j| (i, j)))
        .map(|(i, j)| (matrix[i][j] - matrix[j][i]).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(asym > 1.);
    for (s, before) in old.iter().zip(before) {
        let soil::State::MohrCoulomb(s) = s else {
            panic!()
        };
        assert_eq!(bits(s), before);
    }
    Ok(())
}
#[test]
fn projected_mc_element_all_nine_nonsymmetric_jacobian_blocks_match_independent_residual_differences()
-> Result<()> {
    check_jacobian(Volume::Projected)
}
#[test]
fn independent_volume_mc_element_all_nine_jacobian_blocks_match_independent_residual_differences()
-> Result<()> {
    check_jacobian(Volume::Independent)
}
fn rank(mut a: Vec<Vec<f64>>) -> usize {
    let n = a.len();
    for _ in 0..3 {
        for row in &mut a {
            let scale = row.iter().map(|v| v.abs()).fold(0., f64::max);
            if scale > 0. {
                for v in row {
                    *v /= scale;
                }
            }
        }
        for j in 0..n {
            let scale = a.iter().map(|r| r[j].abs()).fold(0., f64::max);
            if scale > 0. {
                for r in &mut a {
                    r[j] /= scale;
                }
            }
        }
    }
    let mut r = 0;
    for c in 0..n {
        let pivot = (r..n)
            .max_by(|&i, &j| a[i][c].abs().total_cmp(&a[j][c].abs()))
            .unwrap();
        if a[pivot][c].abs() < 1e-9 {
            continue;
        }
        a.swap(r, pivot);
        let d = a[r][c];
        for j in c..n {
            a[r][j] /= d;
        }
        for i in r + 1..n {
            let d = a[i][c];
            for j in c..n {
                a[i][j] -= d * a[r][j];
            }
        }
        r += 1;
        if r == n {
            break;
        }
    }
    r
}
#[test]
fn elastic_both_volume_spaces_have_exactly_six_rigid_modes_and_no_extra_element_pressure_nullspace()
-> Result<()> {
    for volume in [Volume::Projected, Volume::Independent] {
        let e = element(volume, false)?;
        let mut p = parameters(0.);
        p.cohesion_kpa = 1e6;
        let model = soil::Model::MohrCoulomb(p.clone());
        let (old, x) = fixture(&e, &p)?;
        let matrix = e.evaluate(&model, &old, &x, true)?.tangent.unwrap();
        assert_eq!(rank(matrix.clone()), e.count() - 6);
        let scale = matrix.iter().flatten().map(|v| v.abs()).fold(0., f64::max);
        for kind in 0..6 {
            let mut z = vec![0.; e.count()];
            for (i, node) in e.nodes.iter().enumerate() {
                let v = match kind {
                    0 => [1., 0., 0.],
                    1 => [0., 1., 0.],
                    2 => [0., 0., 1.],
                    3 => [0., -node[2], node[1]],
                    4 => [node[2], 0., -node[0]],
                    _ => [-node[1], node[0], 0.],
                };
                z[e.nt + 3 * i..e.nt + 3 * i + 3].copy_from_slice(&v);
            }
            let residual = matrix
                .iter()
                .map(|row| row.iter().zip(&z).map(|(a, b)| a * b).sum::<f64>().abs())
                .fold(0., f64::max);
            assert!(residual / scale < 1e-10);
        }
    }
    Ok(())
}
#[test]
fn curved_element_weak_pressure_and_volume_consistency_does_not_certify_physical_pointwise_stress_equilibrium()
-> Result<()> {
    let e = element(Volume::Projected, true)?;
    let p = parameters(5.);
    let model = soil::Model::MohrCoulomb(p.clone());
    let (old, mut x) = fixture(&e, &p)?;
    x[e.nt + 3 * 4 + 2] += 0.003;
    x[e.nt + 3 * 7] += 0.002;
    let compression: Vec<_> = e
        .points
        .iter()
        .map(|q| {
            -(0..30)
                .map(|j| (q.b[(0, j)] + q.b[(1, j)] + q.b[(2, j)]) * x[e.nt + j])
                .sum::<f64>()
        })
        .collect();
    x[..4].copy_from_slice(&e.project(&compression));
    let first = e.evaluate(&model, &old, &x, false)?;
    x[e.nt + 30..].copy_from_slice(&e.project(&first.mean));
    let final_state = e.evaluate(&model, &old, &x, false)?;
    assert!(final_state.residual[..4].iter().all(|v| v.abs() < 1e-10));
    assert!(
        final_state.residual[e.nt + 30..]
            .iter()
            .all(|v| v.abs() < 1e-10)
    );
    assert!(final_state.maximum_pressure_gap > 1e-5);
    // If a fictitious external force equals the mixed internal force, the
    // mixed u residual vanishes by construction but actual constitutive
    // stresses still give a different physical internal force.
    let gap = final_state
        .actual_internal
        .iter()
        .zip(&final_state.residual[e.nt..e.nt + 30])
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    assert!(
        gap > 1e-6,
        "curved physical-force counterexample not established: {gap}"
    );
    assert!(
        final_state
            .states
            .iter()
            .all(|s| matches!(s,soil::State::MohrCoulomb(v) if p.normalized_yield(v)<=1e-7))
    );
    Ok(())
}
