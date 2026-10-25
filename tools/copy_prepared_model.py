"""Copy a geometry-only snapshot and explicitly replace its physics input."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(source,ground,destination):
    source=source.resolve();manifest=source/'prepared-model.json';data=json.loads(manifest.read_text())
    assert data['schema']=='soil-fem-prepared-model/1' and not data['equilibrium_verified']
    destination.mkdir(parents=True,exist_ok=False)
    for name,digest in data['files_sha256'].items():
        if Path(name).name!=name or '/' in name or '\\' in name:raise ValueError('prepared inventory must contain plain filenames')
        path=(source/name).resolve()
        if path.parent!=source or sha(path)!=digest:raise ValueError(f'prepared file changed or escaped its folder: {name}')
        target=destination/('prepared-ground.json' if name=='ground.json' else name)
        with target.open('xb') as f:f.write(path.read_bytes())
        assert sha(target)==digest
    with (destination/'ground.json').open('xb') as f:f.write(ground.read_bytes())
    record={'schema':'soil-fem-prepared-model-import/1','prepared_source':str(source),'prepared_manifest_sha256':sha(manifest),
        'original_inventory':data,'ground_sha256':sha(destination/'ground.json'),'physics_replaced_explicitly':True}
    (destination/'model-input-package.json').write_text(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('ground',type=Path);p.add_argument('destination',type=Path);a=p.parse_args();copy(a.source,a.ground,a.destination)
