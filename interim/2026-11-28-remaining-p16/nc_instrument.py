"""Select point 351 at the named predictor strain in the two disposable drivers."""
import sys
from pathlib import Path

root = Path(sys.argv[1])


def edit(file, key, new):
    p = root/file
    src = p.read_text()
    assert src.count(key) == 1, (file,key)
    p.write_text(src.replace(key,new))


edit("src/mcc.rs","mod midpoint;",'''mod midpoint;
thread_local! { static P16_POINT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) }; }
pub fn p16_point(index: usize, strain: &[f64;6]) {
    let selected = index==351 && std::env::var("P16_STRAIN").ok().is_some_and(|s| {
        let target: [f64;6] = serde_json::from_str(&s).unwrap();
        strain.iter().zip(target).all(|(x,y)| (*x-y).abs() <= 1e-13*y.abs().max(1e-20))
    });
    P16_POINT.set(selected);
}''')
key = '        let law = &laws[ids[i]];'
edit("examples/p5_homogeneous.rs",key,
     '        soil_fem_3d::mcc::p16_point(i,&strains[i]);\n'+key)
key = '        let reason = law.'
edit("examples/p5_point.rs",key,
     '        soil_fem_3d::mcc::p16_point(row["point"].as_u64().unwrap() as usize,&strain);\n'+key)
key = '        if trace {\n            eprintln!'
p = root/"src/mcc_midpoint.rs"
s = p.read_text()
idx = s.index('        let correction = jac')
pos = s.index('        if trace {',idx)
s = s[:pos]+'''        let correction = if iteration == 0 && super::P16_POINT.get() {
            use std::io::Write;
            use std::process::{Command, Stdio};
            let mut child = Command::new(std::env::var_os("P16_PYTHON").unwrap())
                .arg(std::env::var_os("P16_COMPENSATE").unwrap())
                .stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
            let rows: [[f64;8];8] = std::array::from_fn(|i| std::array::from_fn(|j| jac[(i,j)]));
            let record = serde_json::json!({"matrix":rows,"rhs":value.as_slice(),"native":correction.as_slice()});
            child.stdin.take().unwrap().write_all(serde_json::to_string(&record)?.as_bytes())?;
            let done = child.wait_with_output()?;
            if !done.status.success() { return Err(Error::Material("p16 local direction check".into())); }
            SVector::<f64,8>::from(serde_json::from_slice::<[f64;8]>(&done.stdout)?)
        } else { correction };
'''+s[pos:]
p.write_text(s)
