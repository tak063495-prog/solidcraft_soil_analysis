"""Small real-CLI probes and independent gates for TODO acceptance tests.

No recorded success flags, archived measurements or solver matrices are used
as reference values. Scope refusal is distinct from numerical/input failure.
"""
import json
import copy
import os
from pathlib import Path
import re
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


class CapabilityUnavailable(Exception):
    """Only an exact, documented current feature-scope refusal."""


def read(name):
    return json.loads((ROOT/'examples'/name).read_text(encoding='utf-8'))


def invoke(solver, folder, mesh, ground, case, command='solve-stages'):
    folder.mkdir(parents=True, exist_ok=True)
    paths=[]
    for name, value in [('mesh',mesh), ('ground',ground), ('case',case)]:
        path=folder/(name+'.json')
        path.write_text(json.dumps(value,allow_nan=False),encoding='utf-8')
        paths.append(path)
    output=folder/'run.json'
    assert not output.exists(), 'use a fresh output directory'
    # TimeoutExpired and OS errors intentionally propagate; never xfail/skip.
    result=subprocess.run([str(solver),command,*map(str,paths),str(output)],
                          capture_output=True,text=True,encoding='utf-8',timeout=120)
    (folder/'command.log').write_text(result.stdout+'\n'+result.stderr+f'\nexit_code={result.returncode}\n',encoding='utf-8')
    return result,output,paths


def require_run(result, output, refusal):
    if result.returncode:
        assert result.returncode>0, 'signal termination is not an expected scope refusal'
        assert not output.exists(), 'a refused solve must not publish a Run'
        assert re.search(refusal,result.stderr), result.stderr
        raise CapabilityUnavailable(result.stderr.strip())
    run=json.loads(output.read_text(encoding='utf-8'))
    assert run['complete'] is True, run.get('failure_reason')
    return run


def refused(result,output,pattern):
    assert result.returncode>0, result.stderr
    assert re.search(pattern,result.stderr), result.stderr
    assert not output.exists()


def solver_gates(run):
    """README(a)/(b): accepted nonlinear 1e-8, true linear 1e-10."""
    assert run['complete'] is True
    attempts=run['geostatic']['attempts']+run['attempts']
    accepted=[a for a in attempts if a['accepted']]
    assert accepted, 'an empty accepted history cannot pass'
    linear=[]
    for a in accepted:
        assert np.isfinite(a['relative_force_residual']) and a['relative_force_residual']<=1e-8
        linear.extend(l['true_relative_residual'] for l in a['linear_solves'])
    assert linear and np.isfinite(linear).all() and max(linear)<=1e-10


def contact_water_gates(mesh, ground, case, run):
    """New single-JOINT hydrostatic audit, including declared-tip cases.

    Input geometry/policy independently integrates water loads; saved stress
    reassembly is labelled separately from an independent elastic solution.
    Existing verifier allowlists are not altered or bypassed.
    """
    import fem_reference as fem
    from joint_water_reference import Water
    from joint_ext_reference import constraint_space
    from verify_high_order import prepare, internal, physical
    assert len(case['joint_system']['joints'])==1
    water=Water(mesh,ground,case); contact=water.contact
    definition=case['joint_system']['joints'][0]
    q,offset=prepare(mesh,run['integration_points'])
    for audit in run['transition_audits']:
        stage=next(s for s in case['stages'] if s['name']==audit['stage'])
        if stage['kind']!='joint_birth': continue
        prepared=audit['prepared_state']
        jump=np.einsum('qag,g->qa',contact.c,np.asarray(prepared['displacement_m']).ravel())
        assert prepared['joints']['lifecycle'][0]['generation']>=1
        for i,p in enumerate(prepared['joints']['points']):
            for key in ('plastic_slip_m','accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'):
                assert np.all(np.asarray(p[key])==0.), key
            np.testing.assert_allclose(p['birth']['jump_m'],jump[i],atol=1e-10,rtol=0.)
    for frame in run['frames']:
        s=frame['state'];active=np.asarray(s['active_elements'])
        level=s.get('hydrostatic_water_level_m')
        if level is None: level=ground['initialization']['water_level_m']
        field=water.hydrostatic(active,level)
        np.testing.assert_allclose(s['joint_water']['physical_pressure_kpa'],field['pressure'],atol=1e-8,rtol=0.)
        np.testing.assert_allclose(s['joint_water']['nodal_load_kn'],field['load'].ravel(),atol=1e-10,rtol=0.)
        pore=np.concatenate([water.gamma*np.maximum(level-op['position'][:,2],0.) if water.soil[e] else np.zeros(len(op['weight'])) for e,op in enumerate(q)])
        fi,_=internal(mesh,ground,s['material_states'],active,q,offset,s.get('element_material_ids'),pore)
        force=fi.ravel()
        life=s['joints']['lifecycle'][0]
        if life['active']:
            traction=np.asarray([p['traction_kpa'] for p in s['joints']['points']])
            live=np.asarray([active[p['owner_a']] and active[p['owner_b']] for p in contact.points])
            force+=np.einsum('qag,qa,q->g',contact.c,traction,contact.weight*live)
        load,_=physical(mesh,ground,active,s['surface_pressure_kpa'],case['pressure_on_original_surface_only'],s.get('element_material_ids'),contact)
        load+=s['boundary_factor']*fem.conditions_load(mesh,case.get('boundary',{}),active)
        load+=field['load'];load=load.ravel()
        np.testing.assert_allclose(load,s['applied_nodal_load_kn'],atol=1e-10,rtol=0.)
        basis,_=constraint_space(mesh,case,active,contact,life['tied'])
        assert np.linalg.norm(basis.T@(force-load))/max(np.linalg.norm(load),1.)<1e-10
        assert frame['maximum_constraint_error_m']<=1e-10
        for p in s['joints']['points']:
            if not p.get('birth'): continue
            n=np.asarray(p['normal']);traction=np.asarray(p['traction_kpa']);normal=n@traction
            elastic=normal**2/(2*definition['material']['normal_stiffness_kpa_m'])+np.linalg.norm(traction-normal*n)**2/(2*definition['material']['tangential_stiffness_kpa_m'])
            balance=p['birth']['initial_energy_kpa_m']+p['endpoint_work_kpa_m']-elastic-p['plastic_dissipation_kpa_m']-p['numerical_dissipation_kpa_m']
            assert abs(balance)<=1e-10


def restart_exact(solver,folder,paths,run,stop=1):
    """README(a)/(b): actual separate CLI checkpoint/resume, only clock excluded."""
    checkpoint=folder/'checkpoint.json'; resumed=folder/'resumed.json'
    for args in [('checkpoint-stages',*paths,stop,checkpoint),('resume-stages',*paths,checkpoint,resumed)]:
        p=subprocess.run([str(solver),*map(str,args)],capture_output=True,text=True,encoding='utf-8',timeout=120)
        assert p.returncode==0,p.stderr
    def semantic(v):
        if isinstance(v,dict): return {k:semantic(x) for k,x in v.items() if k!='elapsed_seconds'}
        if isinstance(v,list): return [semantic(x) for x in v]
        return v
    # Canonical JSON preserves signed zeros (ordinary Python equality does not).
    canonical=lambda v:json.dumps(semantic(v),sort_keys=True,allow_nan=False,separators=(',',':'))
    other=json.loads(resumed.read_text(encoding='utf-8'))
    assert other['complete'] is True
    assert canonical(other)==canonical(run)


def rehashed_birth_tampering_refused(solver,folder,paths,run):
    """(a) authenticate actual intermediate birth state after outer rehash.

    Four physical mutations; helper only reseals, product resume must reject.
    A parse/OS/timeout failure of the rehasher is a hard failure, not success.
    """
    helper=solver.parent/'examples'/('joint_ext_checkpoint_rehash.exe' if os.name=='nt' else 'joint_ext_checkpoint_rehash')
    assert helper.is_file(), f'Build first: cargo build --locked --release --examples; missing {helper}'
    count=sum(a['accepted'] for a in run['geostatic']['attempts'])
    birth_names={s['name'] for s in run['case']['stages'] if s['kind']=='joint_birth'}
    found=False
    for a in run['attempts']:
        if not a['accepted']: continue
        count+=1
        if a['stage'] in birth_names:
            found=True;break
    assert found
    original=folder/'birth-checkpoint.json'
    p=subprocess.run([str(solver),'checkpoint-stages',*map(str,paths),str(count),str(original)],capture_output=True,text=True,encoding='utf-8',timeout=120)
    assert p.returncode==0,p.stderr
    saved=json.loads(original.read_text(encoding='utf-8'))
    for mode in ('initial_traction','generation','water_pressure','prepared_stress'):
        changed=copy.deepcopy(saved);state=changed['run']['last_committed']
        if mode=='initial_traction': state['joints']['points'][0]['birth']['initial_traction_kpa'][0]+=1.
        elif mode=='generation': state['joints']['lifecycle'][0]['generation']+=1
        elif mode=='water_pressure': state['joint_water']['physical_pressure_kpa'][0]+=1.
        else: changed['run']['transition_audits'][0]['prepared_state']['material_states'][0]['effective_stress_compression_positive_kpa'][0]+=1.
        bad=folder/(mode+'-bad.json');sealed=folder/(mode+'-sealed.json');output=folder/(mode+'-forbidden.json')
        bad.write_text(json.dumps(changed,allow_nan=False),encoding='utf-8')
        p=subprocess.run([str(helper),str(bad),str(sealed)],capture_output=True,text=True,encoding='utf-8',timeout=120)
        assert p.returncode==0,p.stderr
        rebuilt=json.loads(sealed.read_text(encoding='utf-8'))
        assert rebuilt['run_sha256']!=saved['run_sha256']
        p=subprocess.run([str(solver),'resume-stages',*map(str,paths),str(sealed),str(output)],capture_output=True,text=True,encoding='utf-8',timeout=120)
        (folder/(mode+'-refusal.log')).write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
        assert p.returncode>0,p.stderr
        assert re.search(r'checkpoint|replay|recomput|lifecycle|generation|birth',p.stderr,re.I),p.stderr
        assert not output.exists()
