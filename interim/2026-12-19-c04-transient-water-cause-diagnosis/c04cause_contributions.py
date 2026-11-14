"""Compare global partition-policy contrasts, never substitute a local replay."""
import argparse
from pathlib import Path
import numpy as np
import c04cause_metrics as metrics
from c04cause_observe import partitions
from c04cause_driver import reference,write
import json


def main():
    p=argparse.ArgumentParser();p.add_argument('--controls',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();root=a.controls
    n=metrics.load(root/'native-w32/fields.jsonl');b=metrics.load(root/'adaptive-w32/fields.jsonl');c=metrics.load(root/'boosted-w32/fields.jsonl')
    report=dict(passed=False,global_equilibrium_contrast=True,local_replay_used_as_solution=False,fields={})
    for key in ('u','traction','plastic_slip','accumulated','plastic_dissipation'):
        report['fields'][key]=metrics.budget(n[1][key]-b[1][key],n[1][key]-c[1][key])
    read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    model=reference.Model(read(root/'mesh.json'),read(root/'ground.json'),read(root/'native-w32/case.json'))
    counts=[partitions(root/name,len(model.seeds)) for name in ('native-w32','adaptive-w32','boosted-w32')]
    report['posthoc_partition_vectors_native_vs_boosted_exact']=bool(np.array_equal(counts[0],counts[2]))
    report['posthoc_partition_vectors_native_vs_adaptive_exact']=bool(np.array_equal(counts[0],counts[1]))
    report['partition_vectors_were_not_reference_inputs']=True
    write(a.report,report)

if __name__=='__main__':main()
