"""Bind staged construction, independent solve, restart and legacy evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import meshio
import numpy as np


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run,rust,lint,format_log):
    root=Path(__file__).resolve().parents[1];suite=json.loads((run/'verification.json').read_text());solver=run/'solver/soil-fem-3d.exe';digest=sha(solver)
    assert suite['passed'] and suite['solver_sha256']==digest==sha(root/'target/release/soil-fem-3d.exe')
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust.read_text())));assert count==75
    for p in (rust,lint,format_log):assert p.read_text().rstrip().endswith('Exit code: 0')
    reference=json.loads((run/'independent.json').read_text());assert reference['solver_sha256']==digest and len(reference['frames'])==8
    assert max(f['relative_3d_displacement_error'] for f in reference['frames'])<1e-9
    checkpoint=json.loads((run/'checkpoint/verification.json').read_text());assert checkpoint['solver_sha256']==digest and len(checkpoint['restarts'])==8
    assert all(r['exact_states_frames_material_ids_and_audits'] for r in checkpoint['restarts'])
    files=[solver,rust,lint,format_log,run/'verification.json',run/'independent.json',run/'checkpoint/verification.json']
    for name,expected in checkpoint['sha256'].items():p=run/'checkpoint'/name;assert sha(p)==expected;files.append(p)
    regression=run/'geostatic-regression';geo=json.loads((regression/'verification.json').read_text());assert geo['passed'] and geo['solver_sha256']==digest
    cadlog=regression/'cad-tests.log';assert 'Ran 5 tests' in cadlog.read_text() and cadlog.read_text().rstrip().endswith('Exit code: 0')
    compatible=regression/'mcc-regression.json';old=json.loads(compatible.read_text());assert old['exact_committed_state'] and old['exact_frames']
    files.extend([regression/'verification.json',cadlog,compatible])
    for folder in [run/'cad']+[regression/name for name in ('mc-slope','mcc-slope','free-mcc','slope-pit/analysis','mc-pit','mcc-regression-cad','initial-only')]:
        p=folder/'accepted-stage-run.json';accepted=json.loads(p.read_text());assert accepted['solver_binary_sha256']==digest==sha(folder/'solver/soil-fem-3d.exe')
        for path,expected in accepted['inputs_and_results_sha256'].items():f=Path(path);assert sha(f)==expected;files.append(f)
        files.extend([p,folder/'stage-verification.json'])
    data=json.loads((run/'cad/stage-results.json').read_text());mesh=json.loads((run/'cad/ground/mesh.json').read_text())
    assert data['complete'] and data['schema']=='soil-fem-stage-result/5' and len(data['transition_audits'])==4
    assert len(data['excavation_audits'])==1
    rejected=[a for a in data['attempts'] if not a['accepted']];assert len(rejected)==4 and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    for i,frame in enumerate(data['frames']):
        vtk=run/f'cad/results/stage-{i:03}.vtu';raw=run/f'cad/results/integration-points-{i:03}.vtu'
        v=meshio.read(vtk);q=meshio.read(raw);state=frame['state'];ids=state['element_material_ids'];active=np.asarray(state['active_elements'])
        assert np.array_equal(v.cell_data['material_id'][0],np.asarray(ids)[active])
        assert np.array_equal(q.point_data['material_id'],[ids[p['element_index']] for p in data['integration_points'] if active[p['element_index']]])
        files.extend([vtk,raw])
    # New placement must discard archived coordinates of free replacement nodes.
    assert all(b['placement']=='cad_new_nodes' for b in data['case']['stages'] if b['kind']=='birth')
    assert abs(data['frames'][1]['monitors']['fill_top']['mean_displacement_z_m']-data['frames'][-1]['monitors']['fill_top']['mean_displacement_z_m'])<1e-9
    metrics=[{'stage':f['name'],'active_elements':sum(f['state']['active_elements']),'fill_top':f['monitors']['fill_top'],'relative_force_residual':f['relative_force_residual']} for f in data['frames']]
    files.extend(root/p for p in ('Cargo.lock','Cargo.toml','requirements.lock.txt','sources.lock.json','rust-toolchain.toml','README.md'))
    files.extend(p for name in ('src','tools','scripts','tests','examples','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.ps1','.json','.md'))
    report={'schema':'soil-fem-construction-final-validation/1','passed':True,'rust_tests':count,'cad_tests':5,'solver_sha256':digest,
        'nodes':len(mesh['nodes']),'tetrahedra10':len(mesh['tetrahedra10']),'integration_points':len(data['integration_points']),
        'separate_process_exact_restarts':len(checkpoint['restarts']),'mcc_legacy_exact':True,'metrics':metrics,
        'reference_scope':'Independent full 3D elastic staged stress-history solve; plastic BIRTH/MATSET and MCC/MC mixed behavior covered by Rust theory/restart tests and independent force audits',
        'remaining':'RESET, SRM REAPPLY, zero-pressure MCC/wet BIRTH and other P6 requirements remain incomplete',
        'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(set(files))}}
    with (run/'final-validation.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--rust-log',type=Path,required=True);p.add_argument('--lint-log',type=Path,required=True);p.add_argument('--format-log',type=Path,required=True)
    a=p.parse_args();record(*[v.resolve() for v in (a.run,a.rust_log,a.lint_log,a.format_log)])
