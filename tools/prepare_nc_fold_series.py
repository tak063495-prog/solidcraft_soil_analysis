"""Prepare, but do not solve, the fixed original proportional 1/128,1/256 series."""
import argparse
import hashlib
import json
from pathlib import Path
from probe_nc_midpoint import inputs


def run(output):
    output.mkdir(parents=True,exist_ok=False);records=[]
    for denominator in (128,256):
        folder=output/str(denominator);folder.mkdir();mesh,ground,case=inputs('midpoint')
        case['stages'][2]['initial_increment']=1/denominator
        hashes={}
        for name,value in [('mesh',mesh),('ground',ground),('stages',case)]:
            path=folder/(name+'.json');path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
        records.append(dict(initial_increment=1/denominator,input_sha256=hashes,
                            changed_from_original_128=['stages[2].initial_increment'] if denominator==256 else []))
    (output/'series.json').write_text(json.dumps(dict(solved=False,physical_model_changed=False,profiles=records),indent=2)+'\n',encoding='utf-8')
    return records


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    print(json.dumps(run(p.parse_args().output),indent=2))
