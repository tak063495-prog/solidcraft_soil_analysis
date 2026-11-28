"""Check immutable native geometry/normal/attachment fields from inputs."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
from fixtures import read,write
from reference import Model
from joint_reference import System
from verify_high_order import prepare


def audit(folder):
    m=Model(*(read(folder/(key+'.json')) for key in ('mesh','ground','case')))
    run=read(folder/'native/run.json');assert run['complete']
    System(m.mesh,m.case['joint_system'],layout=run['joint_integration_points'])
    prepare(m.mesh,run['integration_points'])
    maximum_normal=maximum_position=maximum_axis=0.;count=0
    with (folder/'native/accepted.jsonl').open() as states,(folder/'native/structure-reports.jsonl').open() as reports:
        for a,b in itertools.zip_longest(states,reports):
            assert a is not None and b is not None
            state=json.loads(a);report=json.loads(b);assert all(state['active_elements'])
            normal=np.asarray([p['normal'] for p in state['joints']['points']])
            maximum_normal=max(maximum_normal,float(np.max(np.abs(normal-m.system.normal))))
            if m.frame:
                r=report['beams'];maximum_position=max(maximum_position,float(np.max(np.abs(np.asarray(r['node_positions_m'])-m.frame.positions))))
                np.testing.assert_array_equal(r['node_active'],m.frame.active_nodes(state['beams']))
                for raw,(_,_,_,axes,*_) in zip(r['members'],m.frame.data):
                    maximum_axis=max(maximum_axis,float(np.max(np.abs(np.asarray(raw['local_axes'])-axes))))
                for raw,(_,_,kernel) in zip(r.get('plates',[]),m.frame.plate_data):
                    maximum_axis=max(maximum_axis,float(np.max(np.abs(np.asarray(raw['local_axes'])-kernel.axes))))
            if m.tendon:
                for raw,d in zip(report['anchors']['anchors'],m.tendon.data):
                    maximum_position=max(maximum_position,float(np.max(np.abs(np.asarray(raw['point_positions_m'])-d['positions']))))
                    np.testing.assert_allclose(raw['segment_lengths_m'],d['lengths'],atol=1e-12,rtol=1e-12)
            count+=1
    assert maximum_normal<1e-12 and maximum_position<1e-12 and maximum_axis<1e-12
    return dict(accepted_states=count,all_native_volume_and_contact_layouts_verified=True,all_active_domain_unchanged=True,
                maximum_normal_component_difference=maximum_normal,maximum_attachment_position_difference_m=maximum_position,
                maximum_structure_axis_difference=maximum_axis,static_geometry_audit_only=True,counted_as_fresh_solve=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True);p.add_argument('--report',type=Path,required=True);a=p.parse_args()
    result={folder.name:audit(folder) for folder in sorted(a.matrix.iterdir()) if (folder/'comparison.json').is_file()}
    write(a.report,result)


if __name__=='__main__':main()
