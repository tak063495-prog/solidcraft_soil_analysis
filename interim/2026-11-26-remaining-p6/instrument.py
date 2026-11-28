"""Add passive output to the P5 scratch probes. Never edit solver sources in place."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent


def insert(path, old, new):
    source = path.read_text()
    assert source.count(old) == 1, old
    path.write_text(source.replace(old,new))


def run(probe):
    joint = probe/"src/stages_joint_path_error.rs"
    insert(joint,"Some(self.trial(mid, load, settings)?)", "Some(self.trial(mid, load.clone(), settings)?)")
    insert(joint,"            let sample = |s: &Committed|", (HERE/"allocation.rs.inc").read_text()+
           "\n            let sample = |s: &Committed|")
    mcc = probe/"src/mcc_midpoint.rs"
    insert(mcc, '"r":value.as_slice(),', '''
                "jac_rows":(0..8).map(|i|(0..8).map(|j|jac[(i,j)]).collect::<Vec<_>>()).collect::<Vec<_>>(),
                "r":value.as_slice(),''')
    marker = "    if f1 <= 1e-11"
    insert(mcc,marker,'''    if std::env::var_os("P5_POINT_TRACE").is_some() {
        eprintln!("{}",serde_json::json!({"predictor":true,"f0":f0,"f1":f1,"loading":loading,
            "pressure":pressure,"pc":pc,"predictor_stress":trial.stress.map(|s| s.v)}));
    }
'''+marker)


if __name__ == "__main__":
    run(Path(sys.argv[1]))
