"""Recreate the two isolated sources and compare actual formatted build inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from prepare_experiment import prepare

def check(evidence,output):
    output.mkdir(parents=True,exist_ok=False);result={};records=[]
    for label,retain in [('control',False),('retain',True)]:
        dest=output/label;prepare(dest,retain);args=['cargo','fmt','--manifest-path',str(dest/'Cargo.toml')];start=time.monotonic()
        with (output/(label+'-fmt.log')).open('x') as log:r=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
        records.append(dict(command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start));r.check_returncode()
        actual=evidence/'experiment'/label;hashes={}
        for p in sorted(dest.rglob('*')):
            if p.is_file():
                name=p.relative_to(dest).as_posix();data=p.read_bytes()
                if data!=(actual/name).read_bytes():raise ValueError('source reproduction differs: '+label+'/'+name)
                hashes[name]=hashlib.sha256(data).hexdigest()
        result[label]=dict(all_files_exact=True,files_sha256=hashes)
    result['processes']=records
    (output/'source-reproduction.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();check(a.evidence.resolve(),a.output.resolve())
