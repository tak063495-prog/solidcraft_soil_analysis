//! Element matrices for the T5 rigid-space, five-DOF and tying controls.
use soil_fem_3d::{
    Error, Result,
    curved_shell::{Kernel, Tying, geometry::geometries, path::Case},
};
use std::fs::File;
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let [case, mode, output] = args.as_slice() else {
        return Err(Error::Input(
            "CASE tied|raw-membrane|raw-shear NEW_OUTPUT".into(),
        ));
    };
    let case: Case = serde_json::from_reader(File::open(case)?)?;
    let tying = match mode.as_str() {
        "tied" => Tying::default(),
        "raw-membrane" => Tying {
            membrane: false,
            shear: true,
        },
        "raw-shear" => Tying {
            membrane: true,
            shear: false,
        },
        _ => return Err(Error::Input("unknown strain control".into())),
    };
    let mut data = vec![];
    for g in geometries(&case.nodes, &case.faces, &case.directors)? {
        let k = Kernel::new(&g, &case.section, tying)?;
        data.push(
            serde_json::json!({"directors":g.directors.map(<[f64;3]>::from),
            "stiffness_columns":k.stiffness.as_slice(),"drill_columns":k.drill.as_slice()}),
        );
    }
    serde_json::to_writer(File::create_new(output)?, &data)?;
    Ok(())
}
