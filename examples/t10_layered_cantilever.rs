//! SSLS111A. Two 0.4 m offset layers share the six DOFs on the reference plane.
use nalgebra::{DMatrix, DVector, Vector3};
use serde_json::json;
use soil_fem_3d::curved_shell::{Kernel, Q36, Tying, geometry::Geometry, path::solve_supported};
use soil_fem_3d::plates::Section;
use soil_fem_3d::{Error, Result};

// A1 is at the origin; the free edge is x = 10 m for every grid.
struct EdgeElement {
    cut: usize,
    right_of_cut: bool,
    y: [f64; 3],
    at_a1: bool,
    offset: f64,
    dofs: [usize; 36],
    kernel: Kernel,
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let (nx, ny) = match args.as_slice() {
        [] => (10, 5),
        [grid] if grid == "12x6" => (12, 6),
        [grid] if grid == "18x9" => (18, 9),
        _ => return Err(Error::Input("expected no argument, 12x6 or 18x9".into())),
    };
    let diagnostic = !args.is_empty();
    let node = |i, j| j * (2 * nx + 1) + i;
    let count = node(2 * nx, 2 * ny) + 1;
    let mut stiffness = DMatrix::zeros(6 * count, 6 * count);
    let mut load = DVector::zeros(6 * count);
    let mut edges = Vec::new();
    let section = Section {
        id: "SSLS111A layer".into(),
        young_kpa: 200e6,
        poisson: 0.,
        thickness_m: 0.4,
        unit_weight_kn_m3: 0.,
        drilling_factor: 1e-4,
    };
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = (2 * i, 2 * j);
            // Split each rectangular cell. Only one triangle touches A1.
            let triangles = [
                [(0, 0), (2, 0), (0, 2), (1, 0), (1, 1), (0, 1)],
                [(2, 2), (0, 2), (2, 0), (1, 2), (1, 1), (2, 1)],
            ];
            for (side, triangle) in triangles.into_iter().enumerate() {
                let ids = triangle.map(|(a, b)| node(x + a, y + b));
                let points = triangle.map(|(a, b)| {
                    Vector3::new(
                        (x + a) as f64 * 10. / (2 * nx) as f64,
                        (y + b) as f64 * 5. / (2 * ny) as f64,
                        0.,
                    )
                });
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
                    let cut = i + side;
                    if ids[0] == 0 || (diagnostic && cut <= 2) {
                        edges.push(EdgeElement {
                            cut,
                            right_of_cut: side == 0,
                            y: [points[0].y, points[2].y, points[5].y],
                            at_a1: ids[0] == 0,
                            offset,
                            dofs,
                            kernel,
                        });
                    }
                }
            }
        }
        // L12 carries -1 kN/m once over the common quadratic edge.
        for (k, weight) in [1., 4., 1.].into_iter().enumerate() {
            load[6 * node(2 * nx, 2 * j + k) + 2] -= weight * (5. / ny as f64) / 6.;
        }
    }
    let free: Vec<_> = (0..6 * count)
        .filter(|i| (i / 6) % (2 * nx + 1) != 0)
        .collect();
    let reduced = stiffness.select_rows(&free).select_columns(&free);
    let sol = solve_supported(reduced, &load.select_rows(&free))?;
    let mut q = DVector::zeros(6 * count);
    for (i, &dof) in free.iter().enumerate() {
        q[dof] = sol[i];
    }
    let mut layers = Vec::new();
    let mut sections = Vec::new();
    for edge in edges {
        let local = Q36::from_fn(|i, _| q[edge.dofs[i]]);
        let samples = edge.kernel.resultants(&local, edge.offset);
        if edge.at_a1 {
            layers.push(json!({"offset_m":edge.offset,"root_samples":samples}));
        }
        if diagnostic {
            sections.push(json!({"cut":edge.cut,"right_of_cut":edge.right_of_cut,
                "y_m":edge.y,"offset_m":edge.offset,"samples":samples}));
        }
    }
    let mut output = json!({"triangles":4*nx*ny,"tip_dz_m":q[6*node(2*nx,2*ny)+2],"layers":layers});
    if diagnostic {
        let reaction = stiffness * q - load;
        let root_sum = |component: usize| {
            (0..=2 * ny)
                .map(|j| reaction[6 * node(0, j) + component])
                .sum::<f64>()
        };
        output["sections"] = json!(sections);
        output["root_force_z_kn"] = json!(root_sum(2));
        output["root_moment_y_kn_m"] = json!(root_sum(4));
    }
    println!("{output}");
    Ok(())
}
