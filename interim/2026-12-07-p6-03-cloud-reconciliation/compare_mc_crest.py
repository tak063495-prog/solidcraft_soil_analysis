"""One small independent comparison, confined to the original MC column crest."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def compare(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    refpath = HERE/'references/pr24/mcc_reference.py'
    spec = importlib.util.spec_from_file_location('p603_pr24_material_reference', refpath)
    ref = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = ref
    spec.loader.exec_module(ref)
    binary = ROOT/'target/release'/('soil-fem-3d.exe' if os.name == 'nt' else 'soil-fem-3d')
    reportpath = output/'body-plane-strain-ci.json'
    proc = subprocess.run([str(binary), 'plane-strain-2d', 'ci', str(reportpath)], capture_output=True, text=True)
    (output/'body.log').write_text(proc.stdout+proc.stderr, encoding='utf-8')
    if proc.returncode:
        raise RuntimeError(f'body incomplete: exit {proc.returncode}; retain output')
    body = json.loads(reportpath.read_text(encoding='utf-8'))
    level = next(r for s in body['series'] if s['element'] == 'tet10' for r in s['levels'] if r['divisions'] == 4)
    rows = []
    # Existing absolute displacement bound, docs/2D_COMPARISON.ja.md and
    # tools/verify_2d_column_reference.py. Check the crest only, not every DOF.
    for stage, actual in zip(('gravity','load','unload'), level['stages'], strict=True):
        expected = ref.mc_column([.2], stage)['uz'][0]
        error = abs(actual['crest_mean_uz_m'] - expected)
        rows.append(dict(stage=stage,body_crest_m=actual['crest_mean_uz_m'],independent_crest_m=expected,
                         absolute_error_m=error,original_displacement_limit_m=5e-8,operator='<',criterion_met=error<5e-8))
    result = dict(schema='p603-mc-crest-comparison/v1',complete=True,comparison_count=1,
                  scope='dry nonassociated MC .2m original column; Tet10 n4 crest only, three stages',
                  criterion_source='docs/2D_COMPARISON.ja.md:33; tools/verify_2d_column_reference.py bounds/native_u_m',
                  rows=rows,criterion_met=all(r['criterion_met'] for r in rows),
                  source_reference_commit='6ac71b999e718f602a21bfba9ff44af9a787d8c4',
                  source_reference_sha256=hashlib.sha256(refpath.read_bytes()).hexdigest(),
                  binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  body_report_sha256=hashlib.sha256(reportpath.read_bytes()).hexdigest(),
                  tet4_original_failure_preserved=body['column_finest_mesh_five_percent']=='fail',
                  continuum_estimate_is_separate=True,new_all_field_acceptance=False,
                  p6_03_complete=False,p6_complete=False)
    (output/'comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if not result['criterion_met'] or not result['tet4_original_failure_preserved']:
        raise RuntimeError('original criterion not met; result retained')
    return result

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    a=p.parse_args();print(json.dumps(compare(a.output),ensure_ascii=False,indent=2))
