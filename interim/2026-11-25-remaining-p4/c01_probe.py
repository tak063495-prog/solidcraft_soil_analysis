"""Replay archived inputs, retaining the 22 common-anchor fields at N256."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
P2 = HERE.parent / '2026-11-25-remaining5-p2'
P3 = HERE.parent / '2026-11-25-remaining-p3'
ROOT = HERE.parents[1]
sys.path.insert(0, str(P2))
from c01_native import Fields


def run(out, n, alternate):
    out.mkdir(parents=True)
    casefile = (P3/'c01/fixed' if alternate else P2/'c01') / f'case-{n}.json'
    case = json.loads(casefile.read_bytes())
    case['joint_system']['global_step_control']['path_error']['slip_scale_m'] = n**2 * (1e-4/1024**2)
    (out/'input.json').write_text(json.dumps(case))
    (out/'options.json').write_text(json.dumps(dict(full_field_grid=256)))
    cmd = [ROOT/'target/release/examples/c01_methods',ROOT/'examples/joint-water-column-tet10-mesh.json',
           ROOT/'examples/ground-joint-water-mc.json',out/'input.json',out/'result.json','--methods',out/'options.json']
    with Fields(out/'fields',256,'native') as fields, subprocess.Popen(cmd,stdout=subprocess.PIPE) as proc:
        for row in proc.stdout:
            fields.native(*json.loads(row))
        if proc.wait():
            raise RuntimeError(proc.returncode)
        assert fields.count == 1792
    result = json.loads((out/'result.json').read_bytes())
    assert result['complete'] and all(s['attempts']<=2065 for s in result['stages'])


if __name__ == '__main__':
    run(Path(sys.argv[1]),int(sys.argv[2]),sys.argv[3]=='alternate')
