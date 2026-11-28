//! Homogeneous point with the production constant-G update and platen control.
//! Output only points 1..51 at the original output coordinates.
use soil_fem_3d::mcc::{Integration, Parameters, State, invariants};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let initial_bulk = (3.32 - 0.2 * 5_f64.ln()) * 5. / 0.05;
    let model = Parameters {
        critical_state_m: 1.02,
        lambda: 0.2,
        kappa: 0.05,
        poisson: (3. * initial_bulk - 500.) / (2. * (3. * initial_bulk + 250.)),
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.32,
        integration: Integration::BackwardEuler,
    };
    println!("sample,h,load_divisions,p,q,pc,v,axial,lateral,ev,ev_plastic,substeps,plastic_steps");
    for (h, divisions) in [
        (1e-4, 1),
        (5e-5, 1),
        (2.5e-5, 1),
        (1.25e-5, 1),
        (6.25e-6, 1),
        (3.125e-6, 1),
        (3.125e-6, 8),
        (1.5625e-6, 16),
    ] {
        let mut state = State {
            effective_stress_compression_positive_kpa: [5., 5., 5., 0., 0., 0.],
            preconsolidation_pressure_kpa: 40.,
            void_ratio: model.reference_void_ratio(5., 40.)?,
            plastic_strain_compression_positive: [0.; 6],
        };
        let dz = 0.5 / (512 * divisions) as f64;
        let mut lateral = 0.;
        let (mut total, mut plastic) = (0, 0);
        for i in 1..=51 * divisions {
            let (mut lo, mut hi) = (-0.002 / divisions as f64, 0.002 / divisions as f64);
            let mut result = None;
            for _ in 0..60 {
                let x = (lo + hi) / 2.;
                let update =
                    model.integrate_with_shear(&state, [x, x, dz, 0., 0., 0.], h, Some(250.))?;
                let err = update.state.effective_stress_compression_positive_kpa[0] - 5.;
                if err.abs() <= 5e-12 {
                    result = Some((x, update));
                    break;
                }
                if err < 0. {
                    lo = x;
                } else {
                    hi = x;
                }
            }
            let (x, update) = result.ok_or("cell-pressure root did not converge")?;
            lateral += x;
            total += update.substeps;
            plastic += update.plastic_substeps;
            state = update.state;
            if i % divisions == 0 {
                let sample = i / divisions;
                let (p, q) = invariants(&state.effective_stress_compression_positive_kpa);
                let pc = state.preconsolidation_pressure_kpa;
                let v = 1. + state.void_ratio;
                let axial = dz * i as f64;
                let ev = axial + 2. * lateral;
                let ep: f64 = state.plastic_strain_compression_positive[..3].iter().sum();
                println!(
                    "{sample},{h},{divisions},{p},{q},{pc},{v},{axial},{lateral},{ev},{ep},{total},{plastic}"
                );
            }
        }
    }
    Ok(())
}
