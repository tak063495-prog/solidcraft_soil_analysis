//! Read-only replay of homogeneous drained loading against the production constant-G law.
use soil_fem_3d::mcc::{Integration, Parameters, State, invariants};
use std::{
    error::Error,
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    let out = Path::new(&args[1]);
    let mut samples = BufWriter::new(File::create_new(out.join("samples.csv"))?);
    let mut trace = BufWriter::new(File::create_new(out.join("plastic.csv"))?);
    writeln!(
        samples,
        "sample,h,divisions,p,q,pc,v,axial,lateral,ev,ep,steps,plastic,control_steps"
    )?;
    writeln!(
        trace,
        "sample,h,p,q,pc,v_before,v_after,ep,sp,log_pc_ratio,p_before,total_ev"
    )?;
    let bulk = (3.32 - 0.2 * 5_f64.ln()) * 5. / 0.05;
    let model = Parameters {
        critical_state_m: 1.02,
        lambda: 0.2,
        kappa: 0.05,
        poisson: (3. * bulk - 500.) / (2. * (3. * bulk + 250.)),
        reference_pressure_kpa: 1.,
        specific_volume_ncl_at_reference: 3.32,
        integration: Integration::BackwardEuler,
    };
    for (h, divisions, end) in [
        (1e-4, 1, 512),
        (5e-5, 1, 512),
        (2.5e-5, 1, 512),
        (1.25e-5, 1, 512),
        (1e-4, 4, 51),
        (5e-5, 4, 51),
        (2.5e-5, 4, 51),
        (1.25e-5, 4, 51),
    ] {
        let mut state = State {
            effective_stress_compression_positive_kpa: [5., 5., 5., 0., 0., 0.],
            preconsolidation_pressure_kpa: 40.,
            void_ratio: model.reference_void_ratio(5., 40.)?,
            plastic_strain_compression_positive: [0.; 6],
        };
        let dz = 0.5 / (512 * divisions) as f64;
        let mut lateral = 0.;
        let (mut total, mut plastic, mut controls) = (0, 0, 0);
        for i in 1..=end * divisions {
            let (mut lo, mut hi) = (-0.002 / divisions as f64, 0.002 / divisions as f64);
            for (x, side) in [(lo, -1.), (hi, 1.)] {
                let u =
                    model.integrate_with_elastic_split(&state, [x, x, dz, 0., 0., 0.], h, 250.)?;
                controls += u.substeps;
                assert!(side * (u.state.effective_stress_compression_positive_kpa[0] - 5.) >= 0.);
            }
            let mut result = None;
            for _ in 0..60 {
                let x = (lo + hi) / 2.;
                let u =
                    model.integrate_with_elastic_split(&state, [x, x, dz, 0., 0., 0.], h, 250.)?;
                controls += u.substeps;
                let err = u.state.effective_stress_compression_positive_kpa[0] - 5.;
                if err.abs() <= 5e-12 {
                    result = Some((x, u));
                    break;
                }
                if err < 0. { lo = x } else { hi = x }
            }
            let (x, u) = result.ok_or("pressure bracket failed")?;
            if divisions == 1 && i <= 50 && h >= 5e-5 {
                let n = ((2. * x.abs() + dz) / h).ceil() as usize;
                let de = [x, x, dz, 0., 0., 0.].map(|y| y / n as f64);
                let mut parent = state.clone();
                for _ in 0..n {
                    let next = model
                        .integrate_with_elastic_split(&parent, de, 1e-3, 250.)?
                        .state;
                    let dep = std::array::from_fn::<_, 3, _>(|j| {
                        next.plastic_strain_compression_positive[j]
                            - parent.plastic_strain_compression_positive[j]
                    });
                    let ep = dep.iter().sum::<f64>();
                    let sp = dep[2] - ep / 3.;
                    if sp != 0. {
                        let (p, q) = invariants(&next.effective_stress_compression_positive_kpa);
                        let pc = next.preconsolidation_pressure_kpa;
                        writeln!(
                            trace,
                            "{i},{h},{p},{q},{pc},{},{},{ep},{sp},{},{},{}",
                            1. + parent.void_ratio,
                            1. + next.void_ratio,
                            (pc / parent.preconsolidation_pressure_kpa).ln(),
                            invariants(&parent.effective_stress_compression_positive_kpa).0,
                            de[..3].iter().sum::<f64>()
                        )?;
                    }
                    parent = next;
                }
                assert!((parent.void_ratio - u.state.void_ratio).abs() < 1e-12);
                assert!(
                    (parent.preconsolidation_pressure_kpa - u.state.preconsolidation_pressure_kpa)
                        .abs()
                        < 1e-10
                );
            }
            total += u.substeps;
            plastic += u.plastic_substeps;
            state = u.state;
            lateral += x;
            if i % divisions == 0 {
                let (p, q) = invariants(&state.effective_stress_compression_positive_kpa);
                let ax = dz * i as f64;
                writeln!(
                    samples,
                    "{},{h},{divisions},{p},{q},{},{},{ax},{lateral},{},{},{total},{plastic},{controls}",
                    i / divisions,
                    state.preconsolidation_pressure_kpa,
                    1. + state.void_ratio,
                    ax + 2. * lateral,
                    state.plastic_strain_compression_positive[..3]
                        .iter()
                        .sum::<f64>()
                )?;
            }
        }
    }
    Ok(())
}
