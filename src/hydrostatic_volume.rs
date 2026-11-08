//! Exact -integral grad(N) p for affine Tet4/Tet10 hydrostatic water.
//! Material quadrature and displayed GP pressures are deliberately separate.
use crate::{Error, Result, elements};
use nalgebra::{Matrix3, Vector3};

/// None means a curved Tet10: retain its existing isoparametric GP rule.
/// The caller must validate the element before taking this affine shortcut.
pub fn affine_load(points: &[[f64; 3]], level: f64, gamma: f64) -> Result<Option<Vec<[f64; 3]>>> {
    if ![4, 10].contains(&points.len())
        || !points.iter().flatten().all(|x| x.is_finite())
        || !level.is_finite()
        || !gamma.is_finite()
        || gamma < 0.
    {
        return Err(Error::Input(
            "invalid affine hydrostatic volume data".into(),
        ));
    }
    let v: [Vector3<f64>; 4] = std::array::from_fn(|i| Vector3::from(points[i]));
    let jac = Matrix3::from_columns(&[v[1] - v[0], v[2] - v[0], v[3] - v[0]]);
    let inverse = jac
        .try_inverse()
        .filter(|_| jac.determinant() > 0.)
        .ok_or_else(|| Error::Input("positive affine hydrostatic tetrahedron required".into()))?;
    if points.len() == 10 {
        let scale = (1..4).map(|i| (v[i] - v[0]).norm()).fold(0., f64::max);
        for (k, [i, j]) in elements::EDGES.into_iter().enumerate() {
            if (Vector3::from(points[k + 4]) - (v[i] + v[j]) / 2.).norm()
                > 32. * f64::EPSILON * scale
            {
                return Ok(None);
            }
        }
    }
    let wet: Vec<_> = (0..4).filter(|&i| v[i][2] <= level).collect();
    let dry: Vec<_> = (0..4).filter(|&i| v[i][2] > level).collect();
    let cross = |i: usize, j: usize| {
        let mut p = v[i] + (v[j] - v[i]) * ((level - v[i][2]) / (v[j][2] - v[i][2]));
        p[2] = level;
        p
    };
    // A two-wet-vertex wedge is three tetrahedra. For three wet vertices,
    // subtract the dry corner from the full *linear extension* of p.
    // This differs from the independent reference's clipped-face centroid fan.
    let mut cells = Vec::new();
    match wet.len() {
        0 => {}
        1 => cells.push((
            [
                v[wet[0]],
                cross(wet[0], dry[0]),
                cross(wet[0], dry[1]),
                cross(wet[0], dry[2]),
            ],
            1.,
        )),
        2 => {
            let [a, b, c, d] = [
                cross(wet[0], dry[0]),
                cross(wet[0], dry[1]),
                cross(wet[1], dry[0]),
                cross(wet[1], dry[1]),
            ];
            cells.push(([v[wet[0]], a, b, v[wet[1]]], 1.));
            cells.push(([v[wet[1]], a, b, c], 1.));
            cells.push(([v[wet[1]], b, c, d], 1.));
        }
        3 => {
            cells.push((v, 1.));
            cells.push((
                [
                    v[dry[0]],
                    cross(dry[0], wet[0]),
                    cross(dry[0], wet[1]),
                    cross(dry[0], wet[2]),
                ],
                -1.,
            ));
        }
        _ => cells.push((v, 1.)),
    }
    let g: [Vector3<f64>; 4] = std::array::from_fn(|i| {
        inverse.transpose()
            * if i == 0 {
                Vector3::repeat(-1.)
            } else {
                Vector3::ith(i - 1, 1.)
            }
    });
    let mut integral_p = 0.;
    let mut integral_pl = [0.; 4];
    for (cell, sign) in cells {
        let volume = sign
            * Matrix3::from_columns(&[cell[1] - cell[0], cell[2] - cell[0], cell[3] - cell[0]])
                .determinant()
                .abs()
            / 6.;
        let p = cell.map(|x| gamma * (level - x[2]));
        let l = cell.map(|x| {
            let r = inverse * (x - v[0]);
            [1. - r.sum(), r[0], r[1], r[2]]
        });
        let sum_p = p.iter().sum::<f64>();
        integral_p += volume * sum_p / 4.;
        for (i, value) in integral_pl.iter_mut().enumerate() {
            *value += volume / 20.
                * (sum_p * l.iter().map(|a| a[i]).sum::<f64>()
                    + (0..4).map(|j| p[j] * l[j][i]).sum::<f64>());
        }
    }
    let mut force = vec![[0.; 3]; points.len()];
    for i in 0..4 {
        force[i] = (-g[i]
            * if points.len() == 4 {
                integral_p
            } else {
                4. * integral_pl[i] - integral_p
            })
        .into();
    }
    if points.len() == 10 {
        for (k, [i, j]) in elements::EDGES.into_iter().enumerate() {
            force[k + 4] = (-4. * (integral_pl[i] * g[j] + integral_pl[j] * g[i])).into();
        }
    }
    if !force.iter().flatten().all(|x| x.is_finite()) {
        return Err(Error::Solve(
            "non-finite affine hydrostatic volume force".into(),
        ));
    }
    Ok(Some(force))
}
