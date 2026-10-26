"""Versioned native-viewer fields for staged linear frame and flat plate results.

Coordinates/displacements remain relative to CAD, independent of soil display
RESET. Local beam/plate section forces are transported without nodal smoothing.
This extraction does not itself certify the full soil mechanical solution.
"""
import copy
import numpy as np
from beam_reference import System, verify_state

SCHEMA='soil-fem-view-structures/1'


def finite(value,shape,name):
    try:out=np.asarray(value,dtype=float)
    except (TypeError,ValueError) as error:raise ValueError('invalid structural '+name) from error
    if out.shape!=shape or not np.isfinite(out).all() or np.any(abs(out)>1e100):
        raise ValueError('invalid structural '+name)
    return out


def fields(mesh,run,index):
    cfg=run['case'].get('beam_system');frame=run['frames'][index];state=frame['state']
    if cfg is None:
        if state.get('beams') is not None or frame.get('beam_results') is not None:
            raise ValueError('structural state/report without a configuration')
        return None
    if not isinstance(cfg,dict) or cfg.get('schema') not in ('soil-fem-beams/1','soil-fem-structures/1'):raise ValueError('unknown structural configuration')
    if any(len(cfg.get(name,[]))>10000 for name in ('nodes','members','plates')):raise ValueError('structural display inventory exceeds 10000')
    history=state.get('beams');report=frame.get('beam_results')
    if not isinstance(history,dict) or not isinstance(report,dict):raise ValueError('missing structural history/report')
    nodes=cfg['nodes'];n=len(nodes)
    ids=[node['id'] for node in nodes]
    if any(not isinstance(i,str) or not 1<=len(i)<=256 for i in ids) or len(ids)!=len(set(ids)):
        raise ValueError('invalid or duplicate structural node ids')
    displacement=finite(history['node_displacement_m'],(n,3),'nodal displacement')
    rotation=finite(history['node_rotation_rad'],(n,3),'nodal rotation')
    for name,width,inventory in (('reference_global_kinematics',12,cfg['members']),
            ('plate_reference_global_kinematics',18,cfg.get('plates',[]))):
        values=history.get(name,[])
        finite(np.asarray(values).reshape((0,width)) if not inventory and values==[] else values,(len(inventory),width),name)
    for name,inventory in (('active_members',cfg['members']),('active_plates',cfg.get('plates',[]))):
        values=history.get(name,[])
        if len(values)!=len(inventory) or any(type(v) is not bool for v in values):raise ValueError('invalid structural active inventory')
    model=System(mesh,cfg)
    soil_u=finite(state['displacement_m'],(len(mesh['nodes']),3),'soil displacement')
    try:verified=verify_state(model,history,report,soil_u,state['boundary_factor'])
    except (AssertionError,IndexError,KeyError,np.linalg.LinAlgError) as error:
        raise ValueError('structural report/history fails independent member equilibrium') from error
    for key,values in (('node_positions_m',model.positions),('node_displacement_m',displacement),('node_rotation_rad',rotation)):
        actual=finite(report.get(key),(n,3),key)
        if not np.allclose(actual,values,rtol=1e-11,atol=1e-12):raise ValueError('structural report node field differs: '+key)
    if not np.array_equal(report['node_active'],model.active_nodes(history)):raise ValueError('structural node activation differs')
    lookup={name:i for i,name in enumerate(ids)}
    def inventory(name,definitions,active,reference):
        reports=report.get(name,[]);by_id={r['id']:r for r in reports}
        expected={m['id'] for m,live in zip(definitions,active,strict=True) if live}
        if len(by_id)!=len(reports) or set(by_id)!=expected:raise ValueError('duplicate or missing active structural results')
        if len({m['id'] for m in definitions})!=len(definitions):raise ValueError('duplicate structural member ids')
        return [dict(id=m['id'],nodes=[lookup[node] for node in m['nodes']],section_id=m['section_id'],
            active=live,reference_global_kinematics=copy.deepcopy(rest),result=copy.deepcopy(by_id.get(m['id'])))
            for m,live,rest in zip(definitions,active,reference,strict=True)]
    members=inventory('members',cfg['members'],history['active_members'],history['reference_global_kinematics'])
    plates=inventory('plates',cfg.get('plates',[]),history.get('active_plates',[]),history.get('plate_reference_global_kinematics',[]))
    sections={s['id']:s for s in cfg.get('plate_sections',[])}
    for plate in plates:plate['thickness_m']=sections[plate['section_id']]['thickness_m']
    # Each result includes its exact local axes, beam start/mid/end sections or
    # all three DKT samples. Inactive results stay null, never fabricated zeros.
    output=dict(schema=SCHEMA,nodes=dict(ids=ids,positions_m=model.positions.tolist(),soil_node_ids=model.attach,
        displacement_m=displacement.tolist(),rotation_rad=rotation.tolist(),active=verified['node_active'].tolist()),
        members=members,plates=plates,
        global_resultants={name:copy.deepcopy(report[name]) for name in (
            'external_force_kn','external_moment_origin_kn_m','support_force_kn','support_moment_origin_kn_m')},
        semantics=dict(coordinates='CAD metres; displacement from CAD; soil display RESET is not applied',
            beam_section_components=['axial_N','shear_Vy','shear_Vz','torsion_T','moment_My','moment_Mz'],
            beam_section_units=['kN','kN','kN','kN m','kN m','kN m'],
            beam_samples='start, middle, end; exact local axes and end releases retained',
            plate_components=['xx','yy','xy'],plate_membrane='tension-positive kN/m',
            plate_bending='kN m/m; positive with reported local curvature; xy is engineering twist',
            plate_samples='raw DKT barycentric locations; no nodal averaging',
            inactive='archived node kinematics and installation references retained; result null'))
    return output
