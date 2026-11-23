//! SSLS111A. Two 0.4 m offset layers share the six DOFs on the reference plane.
use nalgebra::{DMatrix, DVector, Vector3};
use serde_json::json;
use soil_fem_3d::Result;
use soil_fem_3d::curved_shell::{Kernel, Q36, Tying, geometry::Geometry, path::solve_supported};
use soil_fem_3d::plates::Section;

const NX: usize = 10;
const NY: usize = 5;

fn node(i: usize, j: usize) -> usize {
    j * (2 * NX + 1) + i
}

fn main() -> Result<()> {
    let count = node(2 * NX, 2 * NY) + 1;
    let mut stiffness = DMatrix::zeros(6 * count, 6 * count);
    let mut load = DVector::zeros(6 * count);
    let mut root = Vec::new();
    let section = Section {
        id: "SSLS111A layer".into(),
        young_kpa: 200e6,
        poisson: 0.,
        thickness_m: 0.4,
        unit_weight_kn_m3: 0.,
        drilling_factor: 1e-4,
    };
    for j in 0..NY {
        for i in 0..NX {
            let (x, y) = (2 * i, 2 * j);
            // Split each original 1 m square. Only one triangle touches A1.
            let triangles = [
                [(0, 0), (2, 0), (0, 2), (1, 0), (1, 1), (0, 1)],
                [(2, 2), (0, 2), (2, 0), (1, 2), (1, 1), (2, 1)],
            ];
            for triangle in triangles {
                let ids = triangle.map(|(a, b)| node(x + a, y + b));
                let points = triangle
                    .map(|(a, b)| Vector3::new((x + a) as f64 / 2., (y + b) as f64 / 2., 0.));
                for offset in [-0.2, 0.2] {
                    let geo = Geometry {
                        x: points,
                        directors: [Vector3::z(); 6],
                        e: offset,
                    };
                    let kernel = Kernel::new(&geo, &section, Tying::default())?;
                    let dofs: [usize; 36] = std::array::from_fn(|k| 6 * ids[k / 6] + k % 6);
                    for (a, &row) in dofs.iter().enumerate() {
                        for (b, &column) in dofs.iter().enumerate() {
                            stiffness[(row, column)] += kernel.stiffness[(a, b)];
                        }
                    }
                    if ids[0] == 0 {
                        root.push((offset, dofs, kernel));
                    }
                }
            }
        }
        // L12 carries -1 kN/m once over the common quadratic edge.
        for (k, weight) in [1., 4., 1.].into_iter().enumerate() {
            load[6 * node(2 * NX, 2 * j + k) + 2] -= weight / 6.;
        }
    }
    let free: Vec<_> = (0..6 * count)
        .filter(|i| (i / 6) % (2 * NX + 1) != 0)
        .collect();
    let reduced = stiffness.select_rows(&free).select_columns(&free);
    let sol = solve_supported(reduced, &load.select_rows(&free))?;
    let mut q = DVector::zeros(6 * count);
    for (i, &dof) in free.iter().enumerate() {
        q[dof] = sol[i];
    }
    let layers: Vec<_> = root
        .into_iter()
        .map(|(e, dofs, kernel)| {
            let local = Q36::from_fn(|i, _| q[dofs[i]]);
            json!({"offset_m": e, "root_samples":kernel.resultants(&local,e)})
        })
        .collect();
    println!(
        "{}",
        json!({"triangles":4*NX*NY,"tip_dz_m":q[6*node(2*NX,2*NY)+2],"layers":layers})
    );
    Ok(())
}
