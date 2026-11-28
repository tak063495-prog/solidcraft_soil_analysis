"""Observe accepted contact states before a trial and accepted reload half-steps."""
import sys
from pathlib import Path

root = Path(sys.argv[1])
path = root / "src/stages_joint_path_error.rs"
src = path.read_text()
key = 'coarse.state.stage_index == 6 && std::env::var_os("P17_RETURN").is_some()'
assert src.count(key) == 1
path.write_text(src.replace(key, '(' + key + ' || coarse.state.stage_index == 5 && std::env::var_os("P18_HALVES").is_some())'))

path = root / "src/stages.rs"
src = path.read_text()
# At this point the inherited width rules have already proposed the actual endpoint.
key = '            let target: Vec<_> = if matches!(stage, Stage::ResetDisplay { .. }) {'
assert src.count(key) == 1
state_key = '        let mut last_work:'
assert src.count(state_key) == 1
src = src.replace(state_key, """        let mut p18_latest = None;
        let mut p18_prior: Option<(f64,Vec<nalgebra::Vector3<f64>>)> = None;
""" + state_key)
rule = (Path(__file__).parent/"contact_rule.rs").read_text()
src = src.replace(key,rule + key)
path.write_text(src)
