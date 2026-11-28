"""Supplemental observation of freshly solved fields; NOT a further solve.

Compare all stored water-cache components and pressure source strings. On the
explicit-all-exterior-water case also check the nodal hydrostatic buoyancy
identity W_exterior + W_joint - I_pore = integral(N gamma_w e_z).
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np


def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();out=a.run.resolve()
    root=Path(__file__).resolve().parents[2];sources=read(out/'source-fingerprints.json')
    assert all(sha(root/f)==h for f,h in sources.items())
    result=dict(observation_only_not_new_solve=True,a03_complete=False,script_sha256=sha(Path(__file__)),cases={})
    for name in ('hydrostatic-h1_5','hydrostatic-h2_0'):
        folder=out/name;mesh=read(folder/'mesh.json');ground=read(folder/'ground.json');xyz=np.asarray(mesh['nodes'])
        gamma=ground['initialization']['water_unit_weight_kn_m3'];explicit=ground['initialization']['geostatic'].get('hydrostatic_face_sets')==['exterior']
        weights={m['id']:m['unit_weight_kn_m3'] for m in ground['materials']};assert set(weights.values())=={20.}
        maxima={k:0. for k in ('pressure_kpa','face_force_kn','nodal_load_kn','resultant_force_kn','resultant_moment_kn_m')};rows=[]
        with (folder/'native/accepted.jsonl').open() as n,(folder/'native-fields.jsonl').open() as nf,(folder/'independent-fields.jsonl').open() as rf:
            for ns,nline,rline in itertools.zip_longest(n,nf,rf):
                assert ns is not None and nline is not None and rline is not None
                state,native,reference=map(json.loads,(ns,nline,rline));assert native['row']==reference['row']
                accounting=native['row']['hydraulic_accounting'];assert accounting==dict(model='prescribed_hydrostatic',flux=None,storage=None,birth_water_transfer=None,balance_acceptance=None)
                cache=state['joint_water'];f={k:np.asarray(v) for k,v in native['fields'].items()};r={k:np.asarray(v) for k,v in reference['fields'].items()}
                assert state.get('pore_pressure_kpa') is None
                assert cache['sources']==['continuous']*len(r['joint_water_pressure'])
                differences=dict(pressure_kpa=np.asarray(cache['physical_pressure_kpa'])-r['joint_water_pressure'],
                    face_force_kn=np.asarray(cache['face_nodal_force_a_kn'])-r['joint_face_water_force'],
                    nodal_load_kn=np.asarray(cache['nodal_load_kn']).reshape(-1,3)-r['joint_water_load'],
                    resultant_force_kn=np.asarray(cache['resultant_force_kn'])-r['joint_water_load'].sum(axis=0),
                    resultant_moment_kn_m=np.asarray(cache['resultant_moment_kn_m'])-np.cross(xyz,r['joint_water_load']).sum(axis=0))
                for key,d in differences.items():
                    assert np.isfinite(d).all();maxima[key]=max(maxima[key],float(np.max(np.abs(d),initial=0.)))
                row=dict(stage=native['row']['stage'],progress=native['row']['progress'],sources_exact=True,
                         active_elements=int(f['active_elements'].sum()),exterior_faces=int(f['exterior_faces'].sum()),
                         wet_exterior_faces=int(f['wet_exterior_faces'].sum()),
                         buoyancy_identity_applicable=explicit,native_buoyancy_identity_maximum_absolute_kn=None,
                         independent_buoyancy_identity_maximum_absolute_kn=None)
                if explicit:
                    for label,fields in (('native',f),('independent',r)):
                        actual=fields['exterior_water_load']+fields['joint_water_load']-fields['pore_internal']
                        expected=-gamma/20*fields['gravity_load']
                        error=float(np.max(np.abs(actual-expected),initial=0.));assert error<1e-10
                        row[label+'_buoyancy_identity_maximum_absolute_kn']=error
                rows.append(row)
        assert len(rows)==16
        result['cases'][name]=dict(states=len(rows),all_sources_exact=True,maximum_absolute_differences=maxima,rows=rows)
    assert all(sha(root/f)==h for f,h in sources.items())
    with (out/'water-cache-audit.json').open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
    print(json.dumps(result,ensure_ascii=False))
if __name__=='__main__':main()
