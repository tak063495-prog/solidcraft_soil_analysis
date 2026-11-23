//! Standalone diagnostic. All production shell operators come from the linked library.
use nalgebra::{DMatrix, DVector, SMatrix, SVector, Vector3};
use serde_json::json;
use soil_fem_3d::curved_shell::{Kernel, Q36, Tying, geometry::Geometry, path::solve_supported};
use soil_fem_3d::plates::Section;
use std::error::Error;

const NX: usize = 18;
const NY: usize = 9;
const STRIDE: usize = 2 * NX + 1;
const COUNT: usize = STRIDE * (2 * NY + 1);
const E: f64 = 200e6; // kN/m²
const H: f64 = 0.8;
const LENGTH: f64 = 10.;
const WIDTH: f64 = 5.;

fn point(id: usize) -> Vector3<f64> {
    Vector3::new(
        (id % STRIDE) as f64 * LENGTH / (2 * NX) as f64,
        (id / STRIDE) as f64 * WIDTH / (2 * NY) as f64,
        0.,
    )
}

fn powers(r: f64, s: f64) -> SVector<f64, 6> {
    SVector::from_row_slice(&[1., r, s, r * r, r * s, s * s])
}

// Independent integration of the actual shear operator on a flat affine triangle.
// Only the shear stiffness is exchanged. Retain the production recovery operator.
fn integrate_shear(kernel: &mut Kernel) {
    let sites = &kernel.samples[..6];
    let v = SMatrix::<f64, 6, 6>::from_fn(|i, j| powers(sites[i].p[0], sites[i].p[1])[j]);
    let interpolate = v.transpose().try_inverse().unwrap();
    let points = nalgebra::Vector4::new(
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    );
    let weights = [
        0.34785484513745385,
        0.6521451548625461,
        0.6521451548625461,
        0.34785484513745385,
    ];
    for layer in kernel.samples.chunks_exact(6) {
        let volume: f64 = layer.iter().map(|p| p.weight).sum();
        let modulus = E / 2. * 5. / 6.;
        for p in layer {
            kernel.stiffness -= p.transverse.transpose() * p.transverse * (modulus * p.weight);
        }
        for i in 0..4 {
            let r = (1. + points[i]) / 2.;
            for j in 0..4 {
                let s = (1. - r) * (1. + points[j]) / 2.;
                let l = interpolate * powers(r, s);
                let b = layer
                    .iter()
                    .zip(l.iter())
                    .fold(SMatrix::<f64, 2, 36>::zeros(), |result, (p, weight)| {
                        result + p.transverse * *weight
                    });
                kernel.stiffness += b.transpose()
                    * b
                    * (modulus * volume * weights[i] * weights[j] * (1. - r) / 2.);
            }
        }
    }
}

fn exact_dof(node: usize, component: usize) -> f64 {
    // Reissner–Mindlin strip, nu=0, signed line traction -1 kN/m.
    let x = point(node).x;
    let d = E * H.powi(3) / 12.;
    let s = (5. / 6.) * (E / 2.) * H;
    match component {
        2 => -(LENGTH * x * x / 2. - x.powi(3) / 6.) / d - x / s,
        4 => (LENGTH * x - x * x / 2.) / d,
        _ => 0.,
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mode = std::env::args().nth(1).ok_or("missing axis")?;
    let offsets: &[f64] = match mode.as_str() {
        "single_section" => &[0.],
        "baseline" | "shear_direct" | "membrane_direct" | "gauss_edge" | "full_clamp"
        | "duffy_shear" | "si_units" => &[-0.2, 0.2],
        _ => return Err("unknown axis".into()),
    };
    let unit = if mode == "si_units" { 1000. } else { 1. };
    let section = Section {
        id: "A1 fixed-grid diagnosis".into(),
        young_kpa: E * unit,
        poisson: 0.,
        thickness_m: H / offsets.len() as f64,
        unit_weight_kn_m3: 0.,
        drilling_factor: 1e-4,
    };
    let mut stiffness = DMatrix::zeros(6 * COUNT, 6 * COUNT);
    let mut load = DVector::zeros(6 * COUNT);
    let mut roots = Vec::new();
    for j in 0..NY {
        for i in 0..NX {
            let corner = 2 * j * STRIDE + 2 * i;
            let faces = [
                [0, 2, 2 * STRIDE, 1, STRIDE + 1, STRIDE],
                [
                    2 * STRIDE + 2,
                    2 * STRIDE,
                    2,
                    2 * STRIDE + 1,
                    STRIDE + 1,
                    STRIDE + 2,
                ],
            ];
            for face in faces {
                let ids = face.map(|d| d + corner);
                let dofs: [usize; 36] = std::array::from_fn(|k| ids[k / 6] * 6 + k % 6);
                for &e in offsets {
                    let g = Geometry {
                        x: ids.map(point),
                        directors: [Vector3::z(); 6],
                        e,
                    };
                    let mut kernel = Kernel::new(
                        &g,
                        &section,
                        Tying {
                            membrane: mode != "membrane_direct",
                            shear: mode != "shear_direct",
                        },
                    )?;
                    if mode == "duffy_shear" {
                        integrate_shear(&mut kernel);
                    }
                    for (a, &r) in dofs.iter().enumerate() {
                        for (b, &c) in dofs.iter().enumerate() {
                            stiffness[(r, c)] += kernel.stiffness[(a, b)];
                        }
                    }
                    if ids.contains(&0) {
                        roots.push((g, ids, dofs, kernel));
                    }
                }
            }
        }
        let mut weights = [1. / 6., 4. / 6., 1. / 6.];
        if mode == "gauss_edge" {
            weights = [0.; 3];
            for t in [-1. / 3_f64.sqrt(), 1. / 3_f64.sqrt()] {
                weights[0] += t * (t - 1.) / 4.;
                weights[1] += (1. - t * t) / 2.;
                weights[2] += t * (t + 1.) / 4.;
            }
        }
        for (k, weight) in weights.into_iter().enumerate() {
            load[6 * ((2 * j + k) * STRIDE + 2 * NX) + 2] -= unit * WIDTH / NY as f64 * weight;
        }
    }
    let free: Vec<_> = (0..6 * COUNT).filter(|k| (k / 6) % STRIDE != 0).collect();
    let q = if mode == "full_clamp" {
        // Same fixed DOFs, imposed in the full symmetric system.
        let mut imposed = stiffness.clone();
        for j in 0..=2 * NY {
            for c in 0..6 {
                let k = 6 * j * STRIDE + c;
                imposed.row_mut(k).fill(0.);
                imposed.column_mut(k).fill(0.);
                imposed[(k, k)] = 1.;
            }
        }
        solve_supported(imposed, &load)?
    } else {
        let sol = solve_supported(
            stiffness.select_rows(&free).select_columns(&free),
            &load.select_rows(&free),
        )?;
        let mut result = DVector::zeros(load.len());
        for (i, &k) in free.iter().enumerate() {
            result[k] = sol[i];
        }
        result
    };
    let reaction = &stiffness * &q - &load;
    let root_reaction = |c: usize| {
        (0..=2 * NY)
            .map(|j| reaction[6 * j * STRIDE + c])
            .sum::<f64>()
            / unit
    };
    // Observations only. Neither vector participates in the load or the solve.
    let mut residuals = serde_json::Value::Null;
    if mode == "baseline" {
        let exact = DVector::from_fn(6 * COUNT, |i, _| exact_dof(i / 6, i % 6));
        let projected_residual = &stiffness * exact - &load;
        let component_max = |v: &DVector<f64>| -> [f64; 6] {
            std::array::from_fn(|c| {
                free.iter()
                    .filter(|&&k| k % 6 == c)
                    .map(|&k| v[k].abs())
                    .fold(0.0, f64::max)
            })
        };
        residuals = json!({"solved":component_max(&reaction),
            "analytical_nodes":component_max(&projected_residual)});
    }
    let mut layers = Vec::new();
    for (g, ids, dofs, kernel) in roots {
        let solved = Q36::from_fn(|i, _| q[dofs[i]]);
        let mut layer = json!({"offset":g.e,"solved":kernel.resultants(&solved,g.e)});
        if mode == "baseline" || mode == "shear_direct" {
            let other = Kernel::new(
                &g,
                &section,
                Tying {
                    membrane: true,
                    shear: mode == "shear_direct",
                },
            )?;
            layer["swap_readout"] = json!(other.resultants(&solved, g.e));
        }
        if mode == "baseline" {
            layer["solution_at_nodes"] = json!(solved.as_slice());
            layer["node_x"] = json!(ids.map(|id| point(id).x));
            for (label, translations, rotations) in [
                ("exact_translations", true, false),
                ("exact_rotations", false, true),
                ("exact_both", true, true),
            ] {
                let v = Q36::from_fn(|i, _| {
                    if (i % 6 < 3 && translations) || (i % 6 >= 3 && rotations) {
                        exact_dof(ids[i / 6], i % 6)
                    } else {
                        solved[i]
                    }
                });
                layer[label] = json!(kernel.resultants(&v, g.e));
            }
        }
        layers.push(layer);
    }
    println!(
        "{}",
        json!({"mode":mode,"layers":layers,"force_units_per_kn":unit,"tip_dz_m":q[6*(COUNT-1)+2],
        "free_component_max_residual":residuals,"root_reactions_kn_knm":[root_reaction(2),root_reaction(4)]})
    );
    Ok(())
}
