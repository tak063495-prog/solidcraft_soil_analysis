"""Audit solved mixed manufactured fields, natural mixed energy and refinement.

The common element assessment's energy_error_relative is the *primal
displacement* energy norm. It is retained as such, not called a mixed energy
error. Here pressure is an independent primary unknown and the natural mixed
error uses deviatoric strain and p^2/K; exact pressure is zero for this problem.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from verify_incompressibility import G


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def assess(run,output):
    if output.exists():raise ValueError('refusing to overwrite mixed assessment')
    base=read(run/'assessment.json');assert base['all_cases_solved'] and len(base['cases'])==9
    for name,digest in base['evidence_sha256'].items():assert sha(run/name)==digest,name
    exact=base['exact_strain_squared_norm'];rows=[];files=[run/'assessment.json',Path(__file__),Path(ref.__file__)]
    for original in base['cases']:
        n=original['subdivisions'];nu=original['poisson'];tag=f'tet10-n{n}-nu{str(nu).replace(".","_")}'
        mp=run/f'tet10-n{n}-mesh.json';cp=run/(tag+'-case.json');rp=run/(tag+'-result.json')
        mesh=read(mp);case=read(cp);native=read(rp);files.extend([mp,cp,rp])
        pressure_energy=0.;geometry_count=0;strain_energy=0.;pressure_constraint2=0.
        bulk=case['materials'][0]['young_kpa']/(3*(1-2*nu));layout=native['integration_points']
        eps=np.asarray(native['strain_engineering']);sig=np.asarray(native['stress_tension_positive_kpa']);p=np.asarray(native['pressure_compression_positive_kpa'])
        np.testing.assert_allclose(sig[:,:3].mean(axis=1),-p,rtol=1e-10,atol=1e-10)
        for element in range(len(mesh['tetrahedra'])):
            q=ref.volume(mesh,element);start=geometry_count;stop=start+len(q['weight'])
            for local,weight in enumerate(q['weight']):
                item=layout[start+local]
                assert item['element_index']==element and item['local_index']==local
                np.testing.assert_allclose(item['position_m'],q['position'][local],rtol=1e-12,atol=1e-12)
                assert abs(item['weight_m3']-weight)<1e-12
            e=eps[start:stop];mean=e[:,:3].sum(axis=1);dev=e.copy();dev[:,:3]-=mean[:,None]/3
            pe=float(q['weight']@np.square(p[start:stop])/bulk);pressure_energy+=pe
            strain_energy+=.5*G*float(q['weight']@(np.square(dev)*[2.,2.,2.,1.,1.,1.]).sum(axis=1))+.5*pe
            pressure_constraint2+=float(q['weight']@np.square(mean+p[start:stop]/bulk));geometry_count=stop
        assert geometry_count==len(layout)==len(p)
        assert abs(strain_energy-native['strain_energy_kn_m'])<1e-10*max(1.,abs(strain_energy))
        dev_error2=max(0.,original['strain_relative_l2']**2*exact-2/3*original['divergence_l2']**2)
        row={k:original[k] for k in ('subdivisions','poisson','displacement_relative_l2','strain_relative_l2','displacement_norm_ratio')}
        row.update(natural_mixed_energy_error_relative=float(np.sqrt((G*dev_error2+pressure_energy)/(G*exact))),
                   primal_displacement_energy_error_relative=original['energy_error_relative'],
                   primary_pressure_l2_over_sqrt_bulk=float(np.sqrt(pressure_energy)),
                   pointwise_volume_relation_l2_m_to_3_over_2=float(np.sqrt(pressure_constraint2)),
                   weak_pressure_equation_relative_residual=native['pressure_equation_relative_residual'],
                   all_integration_owners_positions_weights_and_mean_pressures_verified=True)
        rows.append(row)
    gates=[]
    for nu in (.3,.49,.4999):
        group=sorted([r for r in rows if r['poisson']==nu],key=lambda r:r['subdivisions'])
        gate=all(a['displacement_relative_l2']>b['displacement_relative_l2'] for a,b in zip(group,group[1:]))
        gates.append(dict(name=f'displacement_refinement_nu_{nu}',passed=gate,errors=[r['displacement_relative_l2'] for r in group]))
        gates.append(dict(name=f'finest_displacement_below_3_percent_nu_{nu}',passed=group[-1]['displacement_relative_l2']<.03))
    for n in (2,4,6):
        group=[r for r in rows if r['subdivisions']==n];ratio=max(r['displacement_relative_l2'] for r in group)/min(r['displacement_relative_l2'] for r in group)
        gates.append(dict(name=f'poisson_sensitivity_n_{n}',passed=ratio<1.1,maximum_to_minimum_error_ratio=ratio))
    for n in (2,4):
        path=run/f'tet10-n{n}-nu0_4999-independent.json';data=read(path);files.append(path)
        gates.append(dict(name=f'independent_fields_and_pressure_rank_n_{n}',passed=data['passed'] and data['pressure_coupling_nullity']==1,
                          pressure_coupling_rank=data['pressure_coupling_rank'],pressure_inf_sup_positive_beta=data['pressure_inf_sup_positive_beta']))
    record=dict(passed=all(g['passed'] for g in gates),goal_complete=False,p6_complete=False,
        scope='Finite-bulk linear manufactured Tet10/P1 cube only; 2/4/6 subdivisions and three Poisson ratios',
        gates=gates,cases=rows,
        distinction='Primal displacement energy and natural mixed energy are separate metrics; weak volume constraint does not imply pointwise zero divergence',
        limitations='No claim of stability on arbitrary curved/material meshes, exact incompressibility, plasticity, construction or native UI support',
        sha256={str(p):sha(p) for p in set(files)})
    with output.open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(dict(passed=record['passed'],gates=gates)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);a=p.parse_args();assess(a.run.resolve(),a.output.resolve())
