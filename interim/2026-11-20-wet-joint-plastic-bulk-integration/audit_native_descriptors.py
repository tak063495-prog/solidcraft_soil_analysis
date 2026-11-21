"""Read-only native observable provenance audit, not another fresh solve."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from study import read,write,Model,sha


def audit(out,report):
    mesh=read(out/'mesh.json');ground=read(out/'ground.json');gamma=ground['initialization']['water_unit_weight_kn_m3'];base_level=ground['initialization']['water_level_m'];rows=[]
    for folder in sorted(out.glob('*-native')):
        if not (folder/'fields.jsonl').exists():continue
        case=read(folder/'stages.json');model=Model(mesh,ground,case);run=read(folder/'raw/run.json')
        native_position=np.asarray([p['position_m'] for p in run['integration_points']]);native_weight=np.asarray([p['weight_m3'] for p in run['integration_points']])
        position_error=float(np.max(np.abs(native_position-model.positions)));weight_error=float(np.max(np.abs(native_weight-model.weights)))
        np.testing.assert_allclose(native_position,model.positions,rtol=0,atol=1e-12);np.testing.assert_allclose(native_weight,model.weights,rtol=1e-10,atol=1e-12)
        # Existing independent System checks all native face/node/shape/normal/GP descriptors.
        from joint_reference import System
        System(mesh,case['joint_system'],run['joint_integration_points'])
        accumulated=np.zeros((len(model.seeds),6));maximum_pore=0.;maximum_strain=0.;count=0
        with (folder/'raw/accepted.jsonl').open() as states,(folder/'raw/constitutive.jsonl').open() as diagnostics,(folder/'fields.jsonl').open() as exported:
            for line,diag,value in itertools.zip_longest(states,diagnostics,exported):
                assert line is not None and diag is not None and value is not None
                state=json.loads(line);d=json.loads(diag);fields=json.loads(value)['fields']
                # None means the existing fixed hydrostatic ground policy, not an unknown PDE field.
                assert state.get('pore_pressure_kpa') is None
                level=state.get('hydrostatic_water_level_m');level=base_level if level is None else level;assert level==base_level
                pore=gamma*np.maximum(level-native_position[:,2],0.)
                maximum_pore=max(maximum_pore,float(np.max(np.abs(pore-fields['pore']))));assert np.array_equal(pore,np.asarray(fields['pore']))
                accumulated-=np.asarray(d['strain_increment'])
                maximum_strain=max(maximum_strain,float(np.max(np.abs(accumulated-fields['strain']))))
                np.testing.assert_allclose(accumulated,fields['strain'],rtol=1e-10,atol=1e-12)
                count+=1
        rows.append(dict(route=folder.name,accepted_states=count,native_gp_position_maximum_difference_m=position_error,native_gp_weight_maximum_difference_m3=weight_error,
                         native_diagnostic_increment_accumulation_vs_exported_strain_maximum_difference=maximum_strain,native_hydrostatic_observable_vs_exported_pore_maximum_difference_kpa=maximum_pore,
                         pore_source='native GP positions and saved level or ground hydrostatic fallback; no transient GP array in this case'))
    assert len(rows)==6
    write(report,dict(new_solve=False,raw_file_diagnostic_only=True,all_native_geometry_descriptors_checked=True,routes=rows,
          inputs_sha256={name:sha(out/name) for name in ('mesh.json','ground.json','stages.json')}))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args();audit(a.output.resolve(),a.report.resolve())
