//! Read-only diagnostic of the exact shrunk kernel case, plus contained controls.
use solidcraft_kernel::{BoolOp, boolean, cylinder, measure, sphere};
use std::io::{self, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    let mut err = io::BufWriter::new(io::stderr().lock());
    let cylinder_volume = std::f64::consts::PI * 15.0_f64.powi(2) * 20.0;
    let hemisphere = 2.0 * std::f64::consts::PI * 12.5_f64.powi(3) / 3.0;
    for offset in [0.0, 2.499, 2.5] {
        let a = cylinder([25.0, 0.0, 0.0].into(), [0.0, 0.0, 1.0].into(), 15.0, 20.0)?;
        let b = sphere([25.0 + offset, 0.0, 0.0].into(), 12.5)?;
        for op in [BoolOp::Union, BoolOp::Intersect, BoolOp::Cut] {
            let expected = match op {
                BoolOp::Union => cylinder_volume + hemisphere,
                BoolOp::Intersect => hemisphere,
                BoolOp::Cut => cylinder_volume - hemisphere,
            };
            let start = std::time::Instant::now();
            match boolean(&a, &b, op) {
                Ok(body) => {
                    let (volume, invalid) = match body {
                        Some(body) => (measure(&body)?.volume, body.validity().len()),
                        None => (0.0, 0),
                    };
                    if !volume.is_finite() {
                        return Err("nonfinite kernel diagnostic volume".into());
                    }
                    writeln!(out, "{{\"offset_mm\":{offset},\"operation\":\"{op:?}\",\"status\":\"ok\",\"volume_mm3\":{volume},\"expected_volume_mm3\":{expected},\"invalidity_count\":{invalid},\"elapsed_seconds\":{}}}", start.elapsed().as_secs_f64())?;
                }
                Err(e) => {
                    writeln!(out, "{{\"offset_mm\":{offset},\"operation\":\"{op:?}\",\"status\":\"error\",\"expected_volume_mm3\":{expected},\"elapsed_seconds\":{}}}", start.elapsed().as_secs_f64())?;
                    writeln!(err, "offset_mm={offset} operation={op:?}: {e}")?;
                }
            }
            out.flush()?;
        }
    }
    Ok(())
}
