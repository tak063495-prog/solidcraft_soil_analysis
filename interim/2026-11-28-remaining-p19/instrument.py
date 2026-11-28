"""Extend the p18 accepted-state predictor into band 9 in the disposable probe."""
from pathlib import Path
import sys

path = Path(sys.argv[1]) / "src/stages.rs"
src = path.read_text()
changes = {
    "(0.62..0.6875).contains(&progress)": "(0.55..0.6875).contains(&progress)",
    "&& progress >= 0.625": "&& progress >= 0.5625",
    "proposal = proposal.min(progress+offset);":
        "if i != 19 { proposal = proposal.min(progress+offset); }",
    'if std::env::var_os("P18_CONTROL").is_some() {proposal}':
        'if (progress >= 0.625 && std::env::var_os("P18_CONTROL").is_some())'
        ' || (progress < 0.625 && std::env::var_os("P19_CONTROL").is_some()) {proposal}',
}
for old, new in changes.items():
    assert src.count(old) == 1, old
    src = src.replace(old, new)
path.write_text(src)
