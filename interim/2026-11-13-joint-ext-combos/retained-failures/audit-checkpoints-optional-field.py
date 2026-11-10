"""Independently reassemble every saved accepted state, including holding loads.

The corresponding full run must be complete. A partial checkpoint alone never
establishes completion. This auditor has the same bounded dry-tip/static-water
scope as joint_ext_reference, and does not accept Richards or steady flow.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem
from joint_reference import System
from joint_water_reference import Water
from joint_ext_reference import constraint_space
from verify_high_order import prepare,internal,physical


def audit(evidence):
    rows=[]
    for name in ('tip-history-transfer_force','tip-history-stress_free','wet-hydrostatic'):
        folder=evidence/name
        mesh=json.loads((folder/'mesh.json').read_text());ground=json.loads((folder/'ground.json').read_text())
        full=json.loads((folder/'run.json').read_text());assert full['complete']
        case=full['case'];assert not case.get('beam_system') and not case.get('anchor_system')
        assert len(case['joint_system']['joints'])==1
        definition=case['joint_system']['joints'][0];wet=definition.get('water') is not None
        assert (wet and name=='wet-hydrostatic' and definition['water']['kind']=='continuous') or (not wet and definition.get('crack_tip'))
        water=Water(mesh,ground,case) if wet else None
        system=water.contact if wet else System(mesh,case['joint_system'])
        q,offset=prepare(mesh,full['integration_points'])
        for path in sorted(folder.glob('checkpoint-*.json')):
            run=json.loads(path.read_text())['run'];assert run['case']==case and run['provenance']==full['provenance']
            s=run['last_committed'];active=np.asarray(s['active_elements']);assert s.get('pore_pressure_kpa') is None
            level=s.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])
            assert not wet or (water.soil.all() and level>=np.asarray(mesh['nodes'])[:,2].max())
            pore=np.concatenate([water.gamma*(level-v['position'][:,2]) if wet else np.zeros(len(v['weight'])) for v in q])
            force,_=internal(mesh,ground,s['material_states'],active,q,offset,s.get('element_material_ids'),pore)
            force=force.ravel();life=s['joints']['lifecycle'][0]
            if life['active']:
                live=np.asarray([active[p['owner_a']] and active[p['owner_b']] for p in system.points])
                traction=np.asarray([p['traction_kpa'] for p in s['joints']['points']])
                force+=np.einsum('qag,qa,q->g',system.c,traction,system.weight*live)
            load,_=physical(mesh,ground,active,s['surface_pressure_kpa'],case['pressure_on_original_surface_only'],s.get('element_material_ids'),system)
            load+=s['boundary_factor']*fem.conditions_load(mesh,case.get('boundary',{}),active)
            water_error=0.
            if wet:
                field=water.hydrostatic(active,level);load+=field['load']
                water_error=float(np.max(np.abs(field['load'].ravel()-s['joint_water']['nodal_load_kn'])));assert water_error<1e-10
            load=load.ravel()
            release=np.asarray(s['excavation_release_load_kn'])+np.asarray(s['transition_release_load_kn'])
            applied=np.asarray(s['applied_nodal_load_kn'])
            load_error=float(np.max(np.abs(load+release-applied)));assert load_error<1e-10,(name,path.name,load_error)
            t,o=constraint_space(mesh,case,active,system,life['tied'])
            residual=float(np.linalg.norm(t.T@(force-load-release))/max(np.linalg.norm(applied),1.));assert residual<1e-10,(name,path.name,residual)
            used=sorted({int(n) for e in np.flatnonzero(active) for n in fem.cells(mesh)[e]});dofs=(3*np.asarray(used)[:,None]+np.arange(3)).ravel()
            u=np.asarray(s['displacement_m']).ravel()-s['boundary_factor']*o
            constraint=float(np.max(np.abs(u-t@(t.T@u))[dofs]));assert constraint<1e-10
            rows.append(dict(profile=name,checkpoint=path.name,stage_index=s['stage_index'],progress=s['stage_progress'],reassembled_relative_residual=residual,all_nodal_rhs_error_kn=load_error,joint_water_load_error_kn=water_error,constraint_error_m=constraint))
    assert len(rows)==41
    return dict(passed=True,completed_parent_runs_required=True,checkpoints=rows,maximum_reassembled_relative_residual=max(r['reassembled_relative_residual'] for r in rows),maximum_all_nodal_rhs_error_kn=max(r['all_nodal_rhs_error_kn'] for r in rows))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=audit(a.evidence)
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='checkpoints'},indent=2))
