"""Bind the fresh baseline to frozen C-04 fields and the captured diagnostic sources."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
from c04cause_driver import write


def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path)
    p.add_argument('--controls',type=Path,required=True);p.add_argument('--report',type=Path);p.add_argument('--capture',action='store_true')
    a=p.parse_args();records={}
    if a.capture:
        target=a.controls/'source';target.mkdir()
        hashes={}
        for source in sorted(Path(__file__).parent.glob('*.py')):
            shutil.copy2(source,target/source.name);hashes[source.name]=hashlib.sha256(source.read_bytes()).hexdigest()
        write(target/'sha256.json',hashes)
        return
    if a.baseline is None or a.report is None:p.error('--baseline and --report required unless --capture')
    for old,new in (('native-w32','native-w32'),('independent-w32','adaptive-w32')):
        left=[json.loads(x) for x in (a.baseline/old/'fields.jsonl').read_text(encoding='utf-8').splitlines()]
        right=[json.loads(x) for x in (a.controls/new/'fields.jsonl').read_text(encoding='utf-8').splitlines()]
        same=left==right
        cases=json.loads((a.baseline/old/'case.json').read_text(encoding='utf-8'))==json.loads((a.controls/new/'case.json').read_text(encoding='utf-8'))
        if not same or not cases:raise ValueError('original baseline changed: '+new)
        records[new]=dict(accepted_states=len(left),fields_object_exact=same,input_case_object_exact=cases)
    source=a.controls/'source';hashes=json.loads((source/'sha256.json').read_text(encoding='utf-8'))
    execution_sources=('c04cause_reference.py','c04cause_driver.py','c04cause_batch.py')
    for name in execution_sources:
        if hashlib.sha256((source/name).read_bytes()).hexdigest()!=hashes[name]:raise ValueError('captured source changed')
        if hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()!=hashes[name]:raise ValueError('execution source changed since capture')
    write(a.report,dict(passed=False,old_evidence_rewritten=False,comparisons=records,
                       captured_execution_sources_unchanged=True,execution_sources={n:hashes[n] for n in execution_sources}))

if __name__=='__main__':main()
