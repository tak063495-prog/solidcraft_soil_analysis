"""Install observation and the opt-in band-13 predictor in the disposable p16 probe."""
from pathlib import Path
import sys

root = Path(sys.argv[1])


def insert(path, needle, before):
    path = root / path
    source = path.read_text()
    assert source.count(needle) == 1, (path, needle)
    path.write_text(source.replace(needle, before + needle))


stages = "src/stages.rs"
insert(stages, "        let mut last_slip:", '''
        let p17 = p16 && std::env::var_os("P17_WORK").is_some();
        let mut last_work: Option<(f64, f64)> = None;
''')
insert(stages, "            let nominal = if p14 &&", '''
            let nominal = if p17 && (0.8125..0.875).contains(&progress)
                && let Some((width, ratio)) = last_work {
                let h = nominal - progress;
                let predicted = ratio * (h / width).powi(2);
                let to = if predicted > 1. { progress + h / 2. } else { nominal };
                eprintln!("{}",serde_json::json!({"p17_parent":progress,"attempt":count,
                    "last_width":width,"last_work_ratio":ratio,"nominal":nominal,
                    "predicted_ratio":predicted,"to":to}));
                to
            } else { nominal };
''')
insert(stages, "                    if p16 {", '''
                    if p17 {
                        last_work = joint_path_error.as_ref().map(|a| (next - from,
                            a.maximum_dissipation_error_kpa_m / a.dissipation_budget_kpa_m));
                    }
''')

# Emit actual parent, midpoint and fine equilibrium before the accepted trial moves.
# The Python reader joins these candidates to final acceptances by the full span.
insert("src/stages_joint_path_error.rs", "        let trial = if control.resolve_slip_path", '''
        if coarse.state.stage_index == 6 && std::env::var_os("P17_RETURN").is_some()
            && control.resolve_slip_path && !audit.reject {
            let sample = |s: &Committed| {
                let points = &s.joints.as_ref().unwrap().points;
                serde_json::json!({"u":s.displacement_m, "load":s.applied_nodal_load_kn,
                    "traction":points.iter().map(|p| p.traction_kpa).collect::<Vec<_>>(),
                    "jump":points.iter().map(|p| p.jump_m).collect::<Vec<_>>()})
            };
            println!("{}",serde_json::json!({"span":[start,midpoint_progress,fine.state.stage_progress],
                "half":sample(&half.state),"end":sample(&fine.state),
                "normals":if start==0. {Some(a.points.iter().map(|p| p.normal).collect::<Vec<_>>())} else {None},
                "initial":if start==0. {Some(sample(old))} else {None}}));
        }
''')

# Mirror the original material diagnostics into a file for point 351 alone.
# P16_POINT remains the narrower strain match controlling compensated directions.
insert("src/mcc.rs", "pub fn p16_point(", '''
thread_local! { pub(crate) static P17_TRACE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }
pub fn p17_trace(index: usize, strain: &[f64;6], count: usize) {
    let enabled = index == 351 && std::env::var_os("P17_NC_LOG").is_some();
    P17_TRACE.set(enabled);
    if enabled {
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new().create(true).append(true)
            .open(std::env::var_os("P17_NC_LOG").unwrap()).unwrap();
        writeln!(f,"{}",serde_json::json!({"evaluation":351,"strain":strain,"count":count})).unwrap();
    }
}
''')
insert("examples/p5_homogeneous.rs", "        soil_fem_3d::mcc::p16_point", '''
        soil_fem_3d::mcc::p17_trace(i,&strains[i],counts[i]);
''')
path = root/"src/mcc_midpoint.rs"
s = path.read_text()
key = 'std::env::var_os("P5_POINT_TRACE").is_some()'
assert key in s
s = s.replace(key, '('+key+' || super::P17_TRACE.get())')
# Only the existing p5 trace statements in this module use eprintln.
key = 'use '
pos = next(i for i,line in enumerate(s.splitlines(True)) if line.startswith(key))
lines = s.splitlines(True)
lines.insert(pos, '''
macro_rules! eprintln {
    ($($arg:tt)*) => {{
        if super::P17_TRACE.get() {
            use std::io::Write;
            let mut f = std::fs::OpenOptions::new().create(true).append(true)
                .open(std::env::var_os("P17_NC_LOG").unwrap()).unwrap();
            writeln!(f,$($arg)*).unwrap();
        } else { std::eprintln!($($arg)*); }
    }};
}
''')
path.write_text(''.join(lines))
