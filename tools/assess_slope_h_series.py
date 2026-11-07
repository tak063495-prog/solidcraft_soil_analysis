"""Assess user-PC Griffiths runs prepared by prepare_slope_h_series.py.

No original 2D displacement reference is read or generated. Stress is compared
as volume-weighted means on the fixed 200 original quad prisms, not transferred
GP values or the full stress-field norm. Optional exact displacement overlap is
potentially expensive and is deliberately absent from ordinary CI.
"""
import argparse
from pathlib import Path
import numpy as np
from plastic_spatial_convergence import read, write, sha, richardson
from tet10_overlap import Overlay


def contains_input(actual, expected):
    """Serde adds defaults to saved Run inputs; every supplied value must agree."""
    if isinstance(expected,dict):
        return isinstance(actual,dict) and all(
            contains_input(actual[k],v) if k in actual else v in ([],{}) for k,v in expected.items())
    if isinstance(expected,list):
        return isinstance(actual,list) and len(actual)==len(expected) and all(
            contains_input(a,b) for a,b in zip(actual,expected,strict=True))
    return actual==expected


def assess(folder, output, overlap=False):
    folder=Path(folder); manifest=read(folder/'preparation.json')
    if not manifest['prepared'] or manifest['material_reduction_factor'] != 1.:
        raise ValueError('unreduced prepared slope series required')
    if [v['relative_h_xz'] for v in manifest['levels']] != [1.,.5,.25] or len(
            {v['depth_divisions'] for v in manifest['levels']}) != 1:
        raise ValueError('three ratio-two in-plane levels with fixed depth required')
    levels=[]; meshes=[]; runs=[]; files={}
    for record in manifest['levels']:
        sub=folder/record['folder']
        for name,digest in record['files_sha256'].items():
            if sha(sub/name)!=digest:
                raise ValueError('prepared input changed: '+str(sub/name))
        mesh=read(sub/'mesh.json'); result=read(sub/'result.json'); parents=read(sub/'parent-quads.json')
        if not result['complete']:
            raise ValueError('incomplete slope solve is not convergence or an Fs upper bound')
        if any(not contains_input(result[name],read(sub/f'{name}.json')) for name in ('ground','case')):
            raise ValueError('result was solved with different prepared ground or case')
        if runs and (result['ground']!=runs[0]['ground'] or result['case']!=runs[0]['case']
                     or result['provenance']['solver_binary_sha256']!=runs[0]['provenance']['solver_binary_sha256']):
            raise ValueError('physical inputs or solver differ across h series')
        frame=result['frames'][0]
        audit=[frame['relative_force_residual'],frame['maximum_normalized_yield'],
               frame['maximum_absolute_engineering_strain_from_seed']]
        global_force=float(np.linalg.norm(frame['force_balance_error_kn'])/
                           max(np.linalg.norm(frame['applied_force_kn']),np.finfo(float).tiny))
        if (not np.isfinite(audit+[global_force]).all() or audit[0]>1e-8 or audit[1]>1e-7
                or audit[2]>.05 or global_force>1e-8):
            raise ValueError('original slope physical audit gates failed')
        state=frame['state']; xyz=np.asarray(mesh['nodes'])
        if len(state['displacement_m']) != len(xyz) or len(parents)!=len(mesh['tetrahedra']):
            raise ValueError('result or source-prism inventory differs')
        highest=xyz[xyz[:,2]==xyz[:,2].max()]
        position=[float(highest[:,0].max()),0.,float(highest[:,2].max())]
        monitor=mesh['nodes'].index(position)
        grouped=np.zeros((200,6)); volumes=np.zeros(200)
        for point,material in zip(result['integration_points'],state['material_states'],strict=True):
            element=point['element_index']
            if not state['active_elements'][element]:
                raise ValueError('the prepared slope has no excavation stage')
            parent=parents[element]; weight=point['weight_m3']
            grouped[parent]+=weight*np.asarray(material['effective_stress_compression_positive_kpa'])
            volumes[parent]+=weight
        if np.any(volumes<=0) or abs(float(volumes.sum())-220)>1e-8:
            raise ValueError('fixed source-prism coverage failed')
        if levels and not np.allclose(volumes,levels[0]['prism_volumes_m3'],rtol=0,atol=1e-8):
            raise ValueError('source prism volumes differ')
        levels.append(dict(folder=record['folder'],elements=len(mesh['tetrahedra']),
            monitor_position_m=position,monitor_uz_m=state['displacement_m'][monitor][2],
            prism_volumes_m3=volumes.tolist(),prism_mean_stresses_kpa=(grouped/volumes[:,None]).tolist(),
            physical_force_residual=frame['relative_force_residual']))
        meshes.append(mesh); runs.append(result)
        for name in ('mesh.json','ground.json','case.json','parent-quads.json','result.json'):
            files[str(sub/name)]=sha(sub/name)
    pairs=[]
    for i in range(len(levels)-1):
        a,b=levels[i:i+2]; weight=np.asarray(b['prism_volumes_m3'])
        delta=np.asarray(a['prism_mean_stresses_kpa'])-b['prism_mean_stresses_kpa']
        squared=float(weight@np.sum(delta**2*[1,1,1,2,2,2],axis=1))
        pair=dict(source_prism_mean_stress_rms_difference_kpa=float(np.sqrt(squared/weight.sum())))
        if overlap:
            overlay=Overlay(meshes[i],meshes[i+1])
            left,right=[r['frames'][0]['state'] for r in runs[i:i+2]]
            pair['displacement']=overlay.compare(left['displacement_m'],right['displacement_m'],
                left['active_elements'],right['active_elements'])
        pairs.append(pair)
    report=dict(schema='soil-fem-griffiths-in-plane-assessment/1',p6_complete=False,
        general_3d_spatial_convergence_demonstrated=False,griffiths_slope_five_percent_closed=False,
        material_reduction_factor=1.,levels=levels,pairs=pairs,
        monitor=richardson([v['monitor_uz_m'] for v in levels]),files_sha256=files,
        limitations=['No original 661-node displacement data: no 5% verdict.',
            'Stress differences are fixed source-prism averages, not the full stress-field norm.',
            'Fixed depth, fixed increments, F=1 only; Fs and general 3D convergence remain open.'])
    write(output,report)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path); parser.add_argument('output',type=Path)
    parser.add_argument('--displacement-overlap',action='store_true')
    args=parser.parse_args()
    assess(args.folder,args.output,args.displacement_overlap)
