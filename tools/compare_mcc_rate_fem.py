"""Compare completed independent rate FEM and Rust states without generating either.

This development report is not a full acceptance: reference/Rust refinement,
rollback, restart and spatial convergence are separate required evidence.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import fem_reference as ref


def compare(reference,rust,mesh=None):
    if not rust['complete']:raise ValueError('cannot compare an incomplete Rust solve as a completed result')
    if len(reference['frames'])!=len(rust['frames']):raise ValueError('frame count mismatch')
    reports=[]
    layout_checked=False
    if mesh is not None:
        point_index=0;layout=rust['integration_points']
        for element in range(len(mesh['tetrahedra'])):
            q=ref.volume(mesh,element)
            for local,weight in enumerate(q['weight']):
                item=layout[point_index]
                if item['element_index']!=element or item['local_index']!=local:raise ValueError('Rust integration-point ordering differs from independent mesh')
                if not np.allclose(item['position_m'],q['position'][local],rtol=0.,atol=1e-8) or not np.isclose(item['weight_m3'],weight,rtol=1e-10,atol=1e-12):
                    raise ValueError('Rust integration geometry differs from independent quadrature')
                point_index+=1
        if point_index!=len(layout):raise ValueError('unexpected Rust integration-point count')
        layout_checked=True
    def relative(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
    reference_origin=np.asarray(reference.get('reference_displacement_m',np.zeros_like(reference['frames'][0]['displacement_m'])))
    native_origin=np.asarray(rust.get('reference_displacement_m',np.zeros_like(rust['frames'][0]['state']['displacement_m'])))
    previous_reference=None;previous_native=None
    for rf,nf in zip(reference['frames'],rust['frames']):
        if rf['name']!=nf['name']:raise ValueError('stage names differ')
        actual=nf['state'];mask=np.asarray(rf['active_elements'])
        if not np.array_equal(mask,actual['active_elements']):raise ValueError('active set mismatch')
        state=np.asarray(rf['material_state']);nq=state.shape[1]
        gp_mask=np.repeat(mask,nq);all_ref=state.reshape(-1,14)
        ref_fields=dict(stress=all_ref[:,:6],pc=np.exp(all_ref[:,6]),void_ratio=np.exp(all_ref[:,7])-1,
                        plastic_strain=all_ref[:,8:])
        raw=actual['material_states'];native_fields=dict(
            stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in raw]),
            pc=np.asarray([s['preconsolidation_pressure_kpa'] for s in raw]),
            void_ratio=np.asarray([s['void_ratio'] for s in raw]),
            plastic_strain=np.asarray([s['plastic_strain_compression_positive'] for s in raw]))
        values={}
        for name,expected in ref_fields.items():
            observed=native_fields[name]
            if observed.shape!=expected.shape:raise ValueError('integration-point inventory mismatch')
            values[name]=dict(relative_all_points=relative(observed,expected),
                maximum_absolute_all_points=float(np.max(np.abs(observed-expected))),
                relative_active_points=relative(observed[gp_mask],expected[gp_mask]))
        ru=np.asarray(rf['displacement_m']);nu=np.asarray(actual['displacement_m'])
        row=dict(name=rf['name'],displacement_relative_all_nodes=relative(nu,ru),
            displacement_maximum_component_error_m=float(np.max(np.abs(nu-ru))),
            display_displacement_relative_all_nodes=relative(nu-native_origin,ru-reference_origin),
            display_displacement_maximum_component_error_m=float(np.max(np.abs((nu-native_origin)-(ru-reference_origin)))),
            independent_relative_force_residual=rf['relative_free_force_residual'],
            independent_plastic_points=int(np.count_nonzero(np.linalg.norm(ref_fields['plastic_strain'],axis=1)>1e-10)),
            rust_plastic_points=int(np.count_nonzero(np.linalg.norm(native_fields['plastic_strain'],axis=1)>1e-10)),
            fields=values,reference_diagnostics=rf['diagnostics'])
        if previous_reference is not None:
            reference_increment=ru-previous_reference;native_increment=nu-previous_native
            row.update(stage_displacement_relative_all_nodes=relative(native_increment,reference_increment),
                stage_displacement_maximum_component_error_m=float(np.max(np.abs(native_increment-reference_increment))))
        if mesh is not None:
            connectivity=np.asarray(mesh.get('tetrahedra10') or mesh['tetrahedra']);used=np.unique(connectivity[mask]);monitors={}
            for name,key in rust['case'].get('monitor_sets',{}).items():
                ids=np.intersect1d(mesh['node_sets'][key],used)
                if not len(ids):continue
                z=ru[ids,2];nz=nu[ids,2]
                monitors[name]=dict(active_nodes=len(ids),independent_mean_from_seed_m=float(z.mean()),rust_mean_from_seed_m=float(nz.mean()),
                    independent_minimum_from_seed_m=float(z.min()),rust_minimum_from_seed_m=float(nz.min()),
                    independent_maximum_from_seed_m=float(z.max()),rust_maximum_from_seed_m=float(nz.max()))
                monitors[name].update(independent_mean_from_display_reference_m=float((ru-reference_origin)[ids,2].mean()),
                                      rust_mean_from_display_reference_m=float((nu-native_origin)[ids,2].mean()))
                if previous_reference is not None:
                    monitors[name].update(independent_stage_increment_mean_m=float((ru-previous_reference)[ids,2].mean()),
                                          rust_stage_increment_mean_m=float((nu-previous_native)[ids,2].mean()))
            row['monitors']=monitors
            node,component=np.unravel_index(np.argmax(np.abs(nu-ru)),nu.shape)
            row['maximum_displacement_error_location']=dict(node=int(node),component=int(component),position_m=mesh['nodes'][node])
        reports.append(row);previous_reference=ru;previous_native=nu
    return dict(schema='soil-fem-independent-mcc-rate-comparison/1',development_comparison_only=True,
                plastic_3d_accepted=False,integration_layout_checked_independently=layout_checked,frames=reports)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('rust',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--mesh',type=Path);p.add_argument('--brief',action='store_true');a=p.parse_args()
    report=compare(json.loads(a.reference.read_text()),json.loads(a.rust.read_text()),json.loads(a.mesh.read_text()) if a.mesh else None)
    paths=[a.reference,a.rust,Path(__file__)]+([a.mesh] if a.mesh else [])
    report['sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    with a.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps([{k:r[k] for k in ('name','displacement_relative_all_nodes','displacement_maximum_component_error_m','independent_plastic_points','rust_plastic_points')} for r in report['frames']] if a.brief else report,indent=2))
