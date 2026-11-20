"""Re-evaluate the original A-02 input bytes using the current binaries."""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / 'interim/2026-11-25-joint-birth-projection-a02'
sys.path.insert(0, str(ROOT / 'tools'))
import fem_reference as fem
from joint_reference import System
from joint_ext_reference import constraint_space
from face_reference import solve as independent_solve
from verify_high_order import internal, physical, prepare


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def extract(out):
    manifest = read(OLD / 'evidence-manifest.json')
    archive = out / 'inputs.zip'
    with archive.open('xb') as dest:
        for part in manifest['parts']:
            data = (OLD / part['name']).read_bytes()
            assert hashlib.sha256(data).hexdigest() == part['sha256']
            dest.write(data)
    with archive.open('rb') as stream:
        assert hashlib.file_digest(stream, 'sha256').hexdigest() == manifest['archive_sha256']
    with zipfile.ZipFile(archive) as z:
        for path, entry in manifest['entries'].items():
            if path.startswith(('fresh/refusals/', 'fresh/projection/', 'fresh/transfer_force/', 'fresh/stress_free/', 'fresh/stress_free-cycle/')) and (
                    Path(path).name in ('mesh.json', 'ground.json', 'case.json')
                    or path.endswith('-request.json')):
                data = z.read(path)
                assert hashlib.sha256(data).hexdigest() == entry['sha256']
                dst = out / path
                dst.parent.mkdir(exist_ok=True, parents=True)
                dst.write_bytes(data)
    archive.unlink()


def audit_frames(folder, run):
    mesh, ground, case = (read(folder / (n+'.json')) for n in ('mesh','ground','case'))
    system = System(mesh, case['joint_system'])
    ops, offsets = prepare(mesh,run['integration_points'])
    results = []
    for frame in run['frames']:
        s = frame['state']
        active = np.asarray(s['active_elements'])
        force, _ = physical(mesh, ground, active, s['surface_pressure_kpa'],
                            case['pressure_on_original_surface_only'], joint_system=system)
        force += s['boundary_factor'] * fem.conditions_load(mesh, case['boundary'], active)
        soil, _ = internal(mesh,ground,s['material_states'],active,ops,offsets,s.get('element_material_ids'))
        jp = s['joints']['points']
        live = np.array([r['active'] for r in frame['joint_results']['points']])
        contact = np.einsum('qag,qa,q->g',system.c,[p['traction_kpa'] for p in jp],
                            live * system.weight).reshape(-1,3)
        life = s['joints']['lifecycle']
        # Build the tied-node space from the input geometry and recorded event state.
        regions = {l.get('region', {}).get('face_index'): l for l in life}
        tied_points = [p for p in system.points if regions.get(p['face_index'],life[0])['tied']]
        from types import SimpleNamespace
        basis, _ = constraint_space(mesh,case,active,SimpleNamespace(points=tied_points),True)
        residual = basis.T @ (soil+contact-force).ravel()
        tips = case['joint_system']['joints'][0]['crack_tip']['edges'][0]
        results.append(dict(stage=frame['name'], area_m2=float(system.weight[live].sum()),
                            generations=[l['generation'] for l in life], active_gp=int(live.sum()),
                            native_relative_force_residual=frame['relative_force_residual'],
                            max_load_reassembly_error_kn=float(np.max(np.abs(force.ravel()-s['applied_nodal_load_kn']))),
                            max_free_residual_kn=float(np.max(np.abs(residual),initial=0)),
                            max_tip_contact_kn=float(np.max(np.abs(contact[tips]))),
                            regions=[dict(region=l.get('region'),active=l['active'],generation=l['generation'],
                                          frozen=[(a['generation'],a.get('endpoint_work_kn_m')) for a in l['archived']])
                                     for l in life]))
    return results



def stable(value):
    if isinstance(value,dict):
        return {k:stable(v) for k,v in value.items() if k!='elapsed_seconds'}
    if isinstance(value,list):
        return [stable(v) for v in value]
    return value


def compare_path(folder,solver,run,states,rehash):
    files = [folder/(n+'.json') for n in ('mesh','ground','case')]
    mesh,ground,case = map(read,files)
    system = System(mesh,case['joint_system'])
    ops, offsets = prepare(mesh,run['integration_points'])
    # Initial equilibrium is solved directly. Subsequent states use the four input fractions.
    initial_count = round(1/case['settings']['initial_increment'])
    rows = []
    gp_elements = np.array([p['element_index'] for p in run['integration_points']])
    for index,ref in enumerate(states):
        count = initial_count+index
        if index==0 or ref['progress']==1:
            frame = next(f for f in run['frames'] if f['name']==ref['stage'])
            s = frame['state']
        else:
            path = folder/f'checkpoint-{count}.json'
            subprocess.run([solver,'checkpoint-stages',*files,str(count),path],check=True,
                           capture_output=True,timeout=7200)
            s = read(path)['run']['last_committed']
            if ref['progress']==0.5:
                resumed = folder/f'resumed-{count}.json'
                subprocess.run([solver,'resume-stages',*files,path,resumed],check=True,capture_output=True,timeout=7200)
                assert json.dumps(stable(read(resumed)),sort_keys=True)==json.dumps(stable(run),sort_keys=True)
                resumed.unlink()
            if ref['stage']=='partial_rebirth' and ref['progress']==0.5:
                saved = read(path)
                for key in ('owner_b','endpoint_work_kn_m','generation','initial_traction_kpa'):
                    modified = read(path)
                    js = modified['run']['last_committed']['joints']
                    if key=='owner_b':
                        js['lifecycle'][0]['region'][key] += 1
                    elif key=='endpoint_work_kn_m':
                        js['lifecycle'][0]['archived'][0][key] += 1.
                    elif key=='generation':
                        js['lifecycle'][1][key] += 1
                    else:
                        js['points'][-1]['birth'][key][0] += 1.
                    raw,sealed,invalid = (folder/(key+end) for end in ('-raw.json','-sealed.json','-run.json'))
                    raw.write_text(json.dumps(modified),encoding='utf-8')
                    subprocess.run([rehash,raw,sealed],check=True,capture_output=True,timeout=7200)
                    assert read(sealed)['run_sha256']!=saved['run_sha256']
                    proc = subprocess.run([solver,'resume-stages',*files,sealed,invalid],capture_output=True,timeout=7200)
                    assert proc.returncode==1 and not invalid.exists() and b'accepted-path replay mismatch' in proc.stderr
                    (folder/(key+'-refused.log')).write_bytes(proc.stderr)
                    raw.unlink()
                    sealed.unlink()
            path.unlink()
        nodes = np.unique(fem.cells(mesh)[ref['active']])
        stress = np.array([p['effective_stress_compression_positive_kpa'] for p in s['material_states']])
        q = np.array([p['traction_kpa'] for p in s['joints']['points']])
        w = np.array([p['endpoint_work_kpa_m'] for p in s['joints']['points']])
        def relative(a,b):
            return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
        assert list(ref['active'])==s['active_elements']
        life = s['joints']['lifecycle']
        ledger = {l.get('region',{}).get('face_index'):l for l in life}
        gen = [ledger.get(p['face_index'],life[0])['generation'] for p in system.points]
        assert np.array_equal(gen,ref['generation'])
        points=s['joints']['points']
        soil,_ = internal(mesh,ground,s['material_states'],ref['active'],ops,offsets)
        internal_force = soil.ravel()+np.einsum('qag,qa,q->g',system.c,q,system.weight*ref['live'])
        reaction = internal_force-s['applied_nodal_load_kn']
        reaction_diff = reaction-(ref['internal']-ref['load'])
        q_relative = relative(q,ref['q']) if np.any(ref['q']) else float(np.linalg.norm(q))
        rows.append(dict(stage=ref['stage'],progress=ref['progress'],accepted_count=count,
                         restarted=ref['progress']==0.5,
                         rejected_resealed_changes=4 if (ref['stage']=='partial_rebirth' and ref['progress']==0.5) else 0,
                         u_relative=relative(np.asarray(s['displacement_m'])[nodes],ref['u'][nodes]),
                         stress_relative=relative(stress[ref['active'][gp_elements]],ref['stress']),
                         traction_relative=q_relative,
                         max_internal_difference_kn=float(np.max(np.abs(internal_force-ref['internal']))),
                         max_reaction_difference_kn=float(np.max(np.abs(reaction_diff))),
                         max_plastic_dissipation_difference_kpa_m=float(np.max(np.abs([p['plastic_dissipation_kpa_m'] for p in points]-ref['pd']))),
                         max_numerical_dissipation_difference_kpa_m=float(np.max(np.abs([p['numerical_dissipation_kpa_m'] for p in points]-ref['numerical']))),
                         max_plastic_slip_difference_m=float(np.max(np.abs([p['plastic_slip_m'] for p in points]-ref['slip']))),
                         max_reference_jump_difference_m=float(np.max(np.abs([p.get('birth',{}).get('jump_m',[0.]*3) for p in points]-ref['reference']))),
                         traction_max_difference_kpa=float(np.max(np.abs(q-ref['q']))),
                         work_max_difference_kpa_m=float(np.max(np.abs(w-ref['work']))),
                         load_max_difference_kn=float(np.max(np.abs(s['applied_nodal_load_kn']-ref['load'])))))
    return rows


def whole_face_controls(out, history):
    # Import the archived all-owner comparison unchanged, including its gates.
    import importlib.util
    sys.path.insert(0,str(OLD))
    spec = importlib.util.spec_from_file_location('all_owner_study',OLD/'study.py')
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    results = {}
    for label in ('transfer_force','stress_free','stress_free-cycle'):
        folder = out/'fresh'/label
        steps = round(1/read(folder/'case.json')['settings']['initial_increment'])
        legacy.native(folder,history,steps)
        subprocess.run([sys.executable,OLD/'study.py','--output',out,'--steps',str(steps),
                        '--independent',folder],check=True,timeout=7200)
        fields = legacy.compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl')['all_fields']
        results[label] = dict(accepted_states=read(folder/'independent-result.json')['accepted_states'],
                              relative={k:fields[k]['relative_l2'] for k in ('u','stress','traction')},
                              maximum_kn={k:fields[k]['maximum_absolute'] for k in ('internal','applied_load')},
                              category_mismatches=sum(v.get('mismatches',0) for v in fields.values()))
        row = results[label]
        row['existing_gates_passed'] = (max(row['relative'].values())<1e-7
                                        and max(row['maximum_kn'].values())<1e-10
                                        and row['category_mismatches']==0)
    return results

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--bin-dir',type=Path,default=ROOT/'target/release')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    extract(out)
    bin_dir = args.bin_dir.resolve()
    suffix = '.exe' if os.name == 'nt' else ''
    solver = bin_dir / ('soil-fem-3d'+suffix)
    probe = bin_dir / 'examples' / ('joint_a02_projection_probe'+suffix)
    # A checkpoint fingerprints the executable. Keep both executables fixed
    # through all fresh solves and separate-process restarts.
    fixed = out/'bin'
    fixed.mkdir()
    solver = Path(shutil.copy2(solver,fixed))
    probe = Path(shutil.copy2(probe,fixed))
    history = Path(shutil.copy2(bin_dir/'examples'/('joint_a02_history_probe'+suffix),fixed))
    rehash = Path(shutil.copy2(bin_dir/'examples'/('joint_ext_checkpoint_rehash'+suffix),fixed))
    rows = {'binary_sha256':{}}
    for exe in (solver,probe,history,rehash):
        with exe.open('rb') as stream:
            rows['binary_sha256'][exe.name] = hashlib.file_digest(stream,'sha256').hexdigest()
    for label in ('body-across-tip-transfer_force','body-across-tip-stress_free','partial-face-rebirth',
                  'tensile-transfer','outside-circle-transfer','inactive-owner-transfer_force','inactive-owner-stress_free'):
        folder = out / 'fresh/refusals' / label
        result = folder / 'run.json'
        with (folder/'new.log').open('x') as log:
            proc = subprocess.run([solver,'solve-stages',*(folder/(n+'.json') for n in ('mesh','ground','case')),
                                   result],stdout=log,stderr=subprocess.STDOUT,timeout=7200)
        rows[label] = dict(exit_code=proc.returncode,has_run=result.exists())
        if result.exists():
            run = read(result)
            rows[label].update(complete=run['complete'],reason=run.get('failure_reason'),frames=audit_frames(folder,run))
        else:
            rows[label]['diagnostic'] = (folder/'new.log').read_text(encoding='utf-8').strip()
        if label.startswith('body-across-tip-') or label=='partial-face-rebirth':
            states, decision = independent_solve(*(read(folder/(n+'.json')) for n in ('mesh','ground','case')))
            rows[label]['independent'] = decision
            if result.exists():
                rows[label]['state_comparison'] = compare_path(folder,solver,run,states,rehash)
            else:
                assert decision['rejected'] and 'Coulomb cone' in rows[label]['diagnostic']
        if label in ('tensile-transfer','outside-circle-transfer'):
            assert proc.returncode == 1 and not result.exists() and 'Coulomb cone' in rows[label]['diagnostic']
        print(label,proc.returncode,flush=True)
    folder = out/'fresh/projection'
    system = System(read(folder/'mesh.json'),read(folder/'case.json')['joint_system'])
    matrix = np.einsum('qag,q->gqa',system.c,np.sqrt(system.weight)).reshape(system.c.shape[-1],-1)
    for req in sorted(folder.glob('*-request.json')):
        label = req.name.removesuffix('-request.json')
        result = folder/(label+'-new.json')
        with (folder/(label+'-new.log')).open('x') as log:
            proc = subprocess.run([probe,folder/'mesh.json',folder/'case.json',req,result],
                                  stdout=log,stderr=subprocess.STDOUT,timeout=7200)
        data = read(result)
        row = dict(accepted=data['accepted'],diagnostic=data['error'],exit_code=proc.returncode)
        assert proc.returncode == (0 if row['accepted'] else 1)
        if label == 'good':
            assert row['accepted']
            f = np.asarray(data['requested_residual_kn'])
            x,_,rank,sv = np.linalg.lstsq(matrix,f,rcond=1e-12)
            q = np.array([p['birth']['initial_traction_kpa'] for p in data['result']['state']['points']])
            row.update(max_all_nodal_error_kn=data['result']['maximum_force_error_kn'],
                       independent_rank=int(rank),singular_values=sv.tolist(),
                       max_gp_difference_kpa=float(np.max(np.abs(q-x.reshape(-1,3)/np.sqrt(system.weight)[:,None]))))
        if label in ('tension','outside-circle','balanced-unreproducible-tip'):
            assert not row['accepted'] and data['result'] is None
        rows['projection/'+label] = row
    rows['whole_faces'] = whole_face_controls(out,history)
    (out/'results.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
