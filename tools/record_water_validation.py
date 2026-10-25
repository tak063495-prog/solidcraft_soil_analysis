"""Seal water only after checking sources, raw mechanics and independent data."""
import argparse
import hashlib
import json
from pathlib import Path
import re
from freeze_sources import manifest


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run):
    root=Path(__file__).resolve().parents[1];target=run/'final-validation.json'
    assert not target.exists() and read(run/'source-snapshot.json')==manifest(),'source changed after pre-build freeze'
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json');assert suite['passed'] and suite['solver_sha256']==digest
    profiles=('mc-tet4','mc-tet10','mcc-tet10')
    logs=['freeze','format','rust-tests','lint','build','cad-tests','split-tests','wet-faces','column-reference',
        'cad-mc','independent-3d','cad-mcc','cad-mcc-history','cad-mcc-checkpoint','legacy']
    logs += [name+tail for name in profiles for tail in ('','-history','-checkpoint')]
    for name in logs:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=(run/'rust-tests.log').read_text(encoding='utf-8');count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=138
    for name in ('wet_interior_island_is_detected_when_all_six_surface_nodes_are_dry',
        'curved_partial_face_matches_independent_root_and_all_six_nodal_integrals',
        'intermediate_water_head_rebuilds_actual_pressure_and_does_not_interpolate_endpoint_loads',
        'plastic_water_history_and_partial_excavation_resume_exactly_and_reject_tampering'):
        assert name+' ... ok' in rust,name
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8')
    independent={}
    for name in ('wet-faces','column-reference','independent-3d','legacy'):
        d=read(run/name/'verification.json');assert d['passed'] and d['solver_sha256']==digest;independent[name]=d
    assert len(independent['wet-faces']['cases'])==4 and len(independent['column-reference']['cases'])==12
    assert independent['wet-faces']['negative_cases']==['unresolved_budget','nonfinite_area_scale']
    cases={};restarts=0
    for name in (*profiles,'cad-mc','cad-mcc'):
        folder=run/name;data=read(folder/'stage-results.json');audit=read(folder/'stage-verification.json');accepted=read(folder/'accepted-stage-run.json')
        assert data['complete'] and data['schema']=='soil-fem-stage-result/11' and data['provenance']['solver_binary_sha256']==digest
        assert accepted['solver_binary_sha256']==digest
        for path,expected in accepted['inputs_and_results_sha256'].items():assert sha(Path(path))==expected
        assert audit['vtu_roundtrip']=='passed' and audit['raw_integration_point_roundtrip']=='passed'
        assert max(f['independent_force_residual'] for f in audit['frames'])<=data['case']['settings']['relative_force_tolerance']
        assert all(a['committed_state_unchanged_on_rejection'] for a in data['attempts'] if not a['accepted'])
        assert any(a.get('rejected_trial_water_changed') for a in data['attempts'] if not a['accepted'])
        assert all(a['zero_release_relative_force_residual']<1e-8 for a in data['excavation_audits'])
        if name!='cad-mc':
            history=read(run/(name+'-history')/'verification.json');assert history['passed'] and history['solver_sha256']==digest
            assert history['exact_accepted_path_and_all_states']
            if 'mcc' in name:assert history['rejected_plastic_trial_exercised']
            cp=read(run/(name+'-checkpoint')/'verification.json');assert cp['solver_sha256']==digest
            assert len(cp['restarts'])==len(data['case']['stages'])
            assert all(r['exact_states_frames_material_ids_and_audits'] and r['saved_hydrostatic_water_level_m'] is not None for r in cp['restarts'])
            assert any(r['saved_nonzero_excavation_holding'] for r in cp['restarts'])
            for path,expected in cp['sha256'].items():assert sha(run/(name+'-checkpoint')/path)==expected
            restarts+=len(cp['restarts'])
        if name.startswith('cad'):
            assert sha(folder/'solver/soil-fem-3d.exe')==digest
            for frame in audit['frames']:
                if frame['stage'] not in ('initial_k0','preload','unload','raise_water'):
                    native=next(f for f in data['frames'] if f['name']==frame['stage']);h=native['state']['hydrostatic_water_level_m']
                    expected=-16.*10.*max(h-6.,0.);assert abs(frame['cut_water_force_kn'][2]-expected)<1e-8
        cases[name]={'accepted_increments':sum(a['accepted'] for a in data['attempts']),'rejected_increments':sum(not a['accepted'] for a in data['attempts']),
            'maximum_independent_force_residual':max(f['independent_force_residual'] for f in audit['frames']),
            'frames':[{'name':f['name'],'level_m':f['state']['hydrostatic_water_level_m'],'monitors':f['monitors']} for f in data['frames']]}
    assert restarts==21,restarts
    for row in independent['legacy']['cases']:assert row['exact_states_and_frames']
    report={'schema':'soil-fem-water-final-validation/1','passed':True,'solver_sha256':digest,'rust_tests':count,
        'cad_and_split_tests':7,'separate_process_exact_restarts':restarts,'source_freeze_before_build':True,
        'frozen_source_files':len(read(run/'source-snapshot.json')),'cases':cases,
        'independent_wet_faces':independent['wet-faces'],'independent_column_refinement':independent['column-reference'],
        'independent_linear_3d':independent['independent-3d'],'mcc_and_idh_legacy_exact':True,
        'limitations':'Prescribed hydrostatic final drained states; fixed input total unit weights. Material-point volume quadrature is explicit and separately refined. Independent actual-CAD 3D displacement reference is the MC elastic branch. No automatic phreatic/seepage boundary, Darcy, Richards, consolidation time, saturation-dependent weight or wet JOINT. Initial horizontal K0 restrictions remain; P6 is incomplete.',
        'source_snapshot_sha256':sha(run/'source-snapshot.zip'),
        'sha256':{p.relative_to(root).as_posix():sha(p) for p in sorted(run.rglob('*')) if p.is_file()}}
    with target.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','cases','independent_wet_faces','independent_column_refinement','independent_linear_3d')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();record(a.run.resolve())
