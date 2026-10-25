"""Real-result roundtrips and refusal of inconsistent primary-pressure exports."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import meshio
from export_mechanics import export


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(out):
    root=Path(__file__).resolve().parents[1];out.mkdir(exist_ok=False)
    column=root/'runs/mixed-column-verification-20261025-v1'
    curved=root/'runs/mixed-curved-cad-20261025-v1'
    legacy=root/'runs/tet10-boundary-test-20261023-224808-4f38f160/high-order'
    inputs=[Path(__file__),root/'tools/export_mechanics.py',root/'tools/high_order_io.py']
    cases=[('column',column/'column-mesh.json',column/'column-result.json'),
           ('curved',curved/'mesh-q5.json',curved/'result-q5.json'),
           ('legacy',legacy/'mesh-q5.json',legacy/'mechanics-q5.json')]
    inputs.extend(p for _,mesh,result in cases for p in (mesh,result))
    fingerprints={str(p):sha(p) for p in inputs};reports=[]
    for name,mp,rp in cases:
        export(mp,rp,out/name);native=read(rp);mesh=read(mp);record=read(out/name/'export-verification.json')
        assert record['raw_integration_point_roundtrip']=='passed' and record['displacement_roundtrip']=='passed'
        if name!='legacy':
            pressure=meshio.read(out/name/'pressure-dofs.vtu');raw=native['pressure_dofs']
            assert len(pressure.points)==len(raw)
            expected=np.array([p['compression_positive_kpa'] for p in raw])
            assert np.array_equal(pressure.point_data['pressure_compression_positive_kpa'],expected)
            assert np.array_equal(pressure.point_data['original_node_id'],[p['node_index'] for p in raw])
            assert np.array_equal(pressure.point_data['material_id'],[p['material_id'] for p in raw])
            groups={}
            for row in raw:groups.setdefault(row['node_index'],[]).append(row)
            interfaces=[rows for rows in groups.values() if len(rows)>1]
            assert interfaces,'the pressure-jump export must contain shared corners with distinct material DOFs'
            if name=='column':
                assert len(interfaces)==9
                assert all(abs(a['compression_positive_kpa']-b['compression_positive_kpa'])>.1 for a,b in interfaces)
            continuum=meshio.read(out/name/'mechanics.vtu')
            assert 'volume_average_pressure_compression_positive_kpa' in continuum.cell_data
            assert 'pressure_compression_positive_kpa' not in continuum.point_data
            reports.append(dict(name=name,primary_pressure_dofs=len(raw),shared_corners=len(interfaces),material_pressures_kept_separate=True))
        else:
            assert not (out/name/'pressure-dofs.vtu').exists()
            reports.append(dict(name=name,legacy_fields_preserved=True))
    rejected=[];original=read(column/'column-result.json')
    for mode in ('missing','duplicate','gp_mismatch','mean_stress_mismatch','interface_swap'):
        native=json.loads(json.dumps(original))
        if mode=='missing':native['pressure_dofs'].pop()
        elif mode=='duplicate':native['pressure_dofs'].append(native['pressure_dofs'][0])
        elif mode=='gp_mismatch':native['pressure_compression_positive_kpa'][0]+=.1
        elif mode=='mean_stress_mismatch':native['stress_tension_positive_kpa'][0][0]+=.3
        else:
            seen={};pair=None
            for row in native['pressure_dofs']:
                if row['node_index'] in seen:pair=(seen[row['node_index']],row);break
                seen[row['node_index']]=row
            assert pair is not None
            pair[0]['compression_positive_kpa'],pair[1]['compression_positive_kpa']=pair[1]['compression_positive_kpa'],pair[0]['compression_positive_kpa']
        path=out/(mode+'-invalid.json')
        with path.open('x') as f:json.dump(native,f)
        destination=out/(mode+'-must-not-exist')
        try:export(column/'column-mesh.json',path,destination)
        except ValueError as error:
            assert not destination.exists()
            rejected.append(dict(case=mode,reason=str(error),no_export_created=True))
        else:raise AssertionError('inconsistent export was accepted: '+mode)
    assert fingerprints=={str(p):sha(p) for p in inputs},'export input/source changed during checks'
    report=dict(passed=True,goal_complete=False,p6_complete=False,
        scope='VTU/CSV export only; real mixed column/curved results and legacy general elasticity',
        cases=reports,rejected_inputs=rejected,input_and_source_sha256=fingerprints,
        evidence_sha256={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()},
        limitations='No native SolidCraft mixed-result display or additional physical validation claim')
    with (out/'verification.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:report[k] for k in ('passed','cases','rejected_inputs')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();run(a.output.resolve())
