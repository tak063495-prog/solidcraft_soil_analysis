use crate::mesh::geometry;
use crate::{Error, Result};
use nalgebra::{SMatrix, SVector};

pub type BMatrix = SMatrix<f64, 6, 12>;
pub type Tensor = SMatrix<f64, 6, 6>;
pub type Stiffness = SMatrix<f64, 12, 12>;

// Tension-positive stress; engineering strains [xx,yy,zz,gamma_xy,gamma_yz,gamma_zx].
pub fn constitutive(e: f64, nu: f64) -> Result<Tensor> {
    if !e.is_finite() || e <= 0.0 || !nu.is_finite() || !(-1.0..0.5).contains(&nu) || nu == -1.0 {
        return Err(Error::Input("require E>0 and -1<nu<0.5".into()));
    }
    let mu = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mut d = Tensor::zeros();
    for i in 0..3 {
        for j in 0..3 {
            d[(i, j)] = lambda;
        }
        d[(i, i)] += 2.0 * mu;
        d[(i + 3, i + 3)] = mu;
    }
    Ok(d)
}

pub fn element(points: &[[f64; 3]; 4], d: &Tensor) -> Result<(f64, BMatrix, Stiffness)> {
    let (volume, g) = geometry(points)?;
    let mut b = BMatrix::zeros();
    for (i, v) in g.iter().enumerate() {
        let k = 3 * i;
        b[(0, k)] = v.x;
        b[(1, k + 1)] = v.y;
        b[(2, k + 2)] = v.z;
        b[(3, k)] = v.y;
        b[(3, k + 1)] = v.x;
        b[(4, k + 1)] = v.z;
        b[(4, k + 2)] = v.y;
        b[(5, k)] = v.z;
        b[(5, k + 2)] = v.x;
    }
    Ok((volume, b, volume * b.transpose() * d * b))
}

pub fn stress(b: &BMatrix, d: &Tensor, u: &[f64; 12]) -> [f64; 6] {
    let s = d * b * SVector::<f64, 12>::from_row_slice(u);
    std::array::from_fn(|i| s[i])
}

#[cfg(test)]
mod tests {
    use super::*;
    const P: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

    #[test]
    fn affine_patch_and_rigid_body_modes() -> Result<()> {
        let d = constitutive(1000., 0.25)?;
        let (v, b, k) = element(&P, &d)?;
        assert!((v - 1. / 6.).abs() < 1e-15);
        let mut u = SVector::<f64, 12>::zeros();
        for (i, p) in P.iter().enumerate() {
            u[3 * i] = 0.01 * p[0] + 0.02 * p[1];
            u[3 * i + 1] = -0.03 * p[1];
            u[3 * i + 2] = 0.04 * p[2];
        }
        assert!(
            (b * u - SVector::<f64, 6>::from_row_slice(&[0.01, -0.03, 0.04, 0.02, 0., 0.])).norm()
                < 1e-14
        );
        assert!((k - k.transpose()).norm() < 1e-12);
        for axis in 0..3 {
            for i in 0..4 {
                for j in 0..3 {
                    u[3 * i + j] = f64::from(j == axis);
                }
            }
            assert!((k * u).norm() < 1e-12);
            for (i, p) in P.iter().enumerate() {
                let a = nalgebra::Vector3::from_row_slice(p);
                let mut w = nalgebra::Vector3::zeros();
                w[axis] = 1.;
                let r = w.cross(&a);
                for j in 0..3 {
                    u[3 * i + j] = r[j];
                }
            }
            assert!((k * u).norm() < 1e-12);
        }
        Ok(())
    }

    #[test]
    fn rejects_invalid_geometry_and_material() {
        let mut p = P;
        p.swap(1, 2);
        assert!(geometry(&p).is_err());
        assert!(constitutive(1., 0.5).is_err());
        assert!(constitutive(f64::NAN, 0.3).is_err());
        assert!(constitutive(1., -1.).is_err());
    }
}
