"""Audit actual Cartesian vectors; distinguish reused historical label strings."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run(source,manifest):
    entries=json.loads(manifest.read_text(encoding='utf-8'))['entries']
    records={}
    for folder in sorted(source.iterdir()):
        if not folder.is_dir() or not (folder/'model.json').exists():continue
        controls=json.loads((folder/'controls.json').read_text(encoding='utf-8'))
        binding=controls.get('binding',controls)
        order=binding['order']
        if order not in ['x-first','x-last','y-first','y-last','z-first','z-last']:continue
        for name in ['model.json','controls.json']:
            if sha(folder/name)!=entries[f'study/{folder.name}/{name}']['sha256']:
                raise ValueError('sealed input differs: '+folder.name+'/'+name)
        model=json.loads((folder/'model.json').read_text(encoding='utf-8'))
        holding=np.asarray(model['holding_kn'])
        axis,sequence=order.split('-')
        mask=np.arange(len(holding))%3=='xyz'.index(axis)
        if sequence=='last':mask=~mask
        first=np.where(mask,holding,0.)
        if not np.array_equal(first,np.asarray(model['first_group_kn'])):
            raise ValueError('actual Cartesian first group differs')
        if not np.array_equal(holding-first,np.asarray(model['second_group_kn'])):
            raise ValueError('actual Cartesian second group differs')
        records[folder.name]=dict(actual_order=order,global_nodal_axis=axis,axis_released=sequence,
            first_group_identity_exact=True,second_group_identity_exact=True,
            legacy_force_decomposition_label=binding.get('force_decomposition'),
            legacy_grouping_label=controls.get('grouping'),
            legacy_fixed_labels_do_not_describe_cartesian_vectors='force_decomposition' in binding,
            model_sha256=sha(folder/'model.json'),controls_sha256=sha(folder/'controls.json'))
    if not records:raise ValueError('no Cartesian evidence found')
    return dict(numerical_acceptance=False,checkpoint_payloads_or_bindings_modified=False,
        archive_manifest_sha256=sha(manifest),audit_source_sha256=sha(Path(__file__)),records=records)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','manifest','output']:p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(run(args.source,args.manifest),indent=2)+'\n',encoding='utf-8')
