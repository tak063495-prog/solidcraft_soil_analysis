"""Expand the exact verification inputs into a new directory, without overwrites."""
from pathlib import Path
import argparse,hashlib,json,zipfile

def unpack(output):
    source=Path(__file__).resolve().parent
    manifest=json.loads((source/'input-manifest.json').read_text(encoding='utf-8'))
    expected={'mesh.json','ground.json','mesh-identity.json','stages-coarse.json','stages-fine.json'}
    assert set(manifest)==expected
    with zipfile.ZipFile(source/'inputs.zip') as archive:
        assert set(archive.namelist())==expected and len(archive.namelist())==len(expected)
        contents={name:archive.read(name) for name in expected}
    for name,data in contents.items():
        assert len(data)==manifest[name]['bytes']
        assert hashlib.sha256(data).hexdigest()==manifest[name]['sha256']
    output.mkdir(parents=True,exist_ok=False)
    for label in ('coarse','fine'):
        folder=output/label;(folder/'ground').mkdir(parents=True)
        files={'ground/mesh.json':'mesh.json','ground/ground.json':'ground.json',
               'mesh-identity.json':'mesh-identity.json','stages.json':'stages-'+label+'.json'}
        for destination,name in files.items():
            with (folder/destination).open('xb') as out:out.write(contents[name])
    print('Exact verification inputs created:',output)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path)
    unpack(parser.parse_args().output.resolve())
