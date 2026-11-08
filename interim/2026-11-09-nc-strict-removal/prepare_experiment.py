"""Isolated midpoint repartition experiment from the verified PR #13 head.

No production source, checkpoint, law, force gate or Newton limit is modified.
The control uses the prior independently checked compensated LU. The variant
retains the current evaluated displacement when the material partition grows.
"""
import argparse
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = '6429e8366ad7b02ab546cc9a9f7ff321b693962b'
PRIOR = ROOT/'interim/2026-11-08-nc-compensated-linear-replay'


def prepare(output, retain):
    spec = importlib.util.spec_from_file_location('prior_nc_prepare', PRIOR/'prepare.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.BASE = BASE
    hashes = module.prepare(output, True)
    if retain:
        path = output/'src/stages.rs'
        text = path.read_text()
        old = '                    u = constrained_start.clone();'
        assert text.count(old) == 1
        text = text.replace(old, '''                    // Experimental only: reevaluate the CURRENT iterate using
                    // the larger fixed partition at the next Newton iteration.
                    // Constitutive history still starts from immutable old.
                    // Do not accept the unevaluated candidate or reuse its residual.''')
        old = '        let constrained_start = u.clone();\n'
        assert text.count(old) == 1
        path.write_text(text.replace(old, ''))
    return hashes


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--retain', action='store_true')
    args = parser.parse_args()
    prepare(args.output.resolve(), args.retain)
