"""Report layout-sensitive residuals without reclassifying recorded trials."""
import argparse
from collections import deque
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from linear_evidence_codec import read_rows


def rows(path):
    with gzip.open(path,'rt') as stream:
        for line in stream:yield json.loads(line)


def audit(evidence):
    result=dict(passed=False,existing_linear_gate=1e-10,all_actual_matrices_available=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),paths={})
    for folder in sorted(evidence.iterdir()):
        if not (folder/'report.json').exists() or not (folder/'controls.json').exists():continue
        systems=iter(rows(folder/'linear-systems.jsonl.gz'));tail=deque(maxlen=128)
        total=0;accepted=0;canonical_over=0;original_over=0;maximum_canonical=0.;maximum_original=0.
        for attempt in rows(folder/'attempts.jsonl.gz'):
            events=attempt['events'];rate=events and events[0].get('method','').startswith('global_load_rk4')
            count=events[0]['actual_linear_systems'] if rate else len(events)
            label=attempt.get('label','step-'+str(attempt.get('step')))
            for _ in range(count):
                row=next(systems)
                if row['label']!=label:raise ValueError('system/attempt alignment differs')
                total+=1;tail.append((row,attempt['accepted']))
                if attempt['accepted']:
                    accepted+=1;a=row['true_relative_residual'];b=row['independently_recomputed_true_relative_residual']
                    maximum_original=max(maximum_original,a);maximum_canonical=max(maximum_canonical,b)
                    original_over+=a>1e-10;canonical_over+=b>1e-10
        if next(systems,None) is not None:raise ValueError('unassigned system summaries')
        retained=list(read_rows(folder/'terminal-linear-systems.jsonl.f64.gz'))
        if len(retained)!=len(tail):raise ValueError('retained system count differs')
        extended_maximum=0.;extended_over=0;retained_accepted=0;f_matches=0;examples=[]
        for actual,(summary,is_accepted) in zip(retained,tail):
            raw=(json.dumps(actual,allow_nan=False)+'\n').encode()
            if hashlib.sha256(raw).hexdigest()!=summary['raw_record_sha256']:raise ValueError('retained system hash differs')
            A=np.array(actual['matrix']);b=np.array(actual['rhs']);x=np.array(actual['solution'])
            fortran=float(np.linalg.norm(np.asfortranarray(A)@x-b)/max(np.linalg.norm(b),1e-30))
            f_matches+=fortran==summary['true_relative_residual']
            al=A.astype(np.longdouble);bl=b.astype(np.longdouble);xl=x.astype(np.longdouble)
            residual=al@xl-bl
            extended=float(np.sqrt(np.sum(residual*residual))/max(np.sqrt(np.sum(bl*bl)),np.longdouble(1e-30)))
            if is_accepted:
                retained_accepted+=1;extended_maximum=max(extended_maximum,extended);extended_over+=extended>1e-10
                if summary['independently_recomputed_true_relative_residual']>1e-10:
                    examples.append(dict(label=summary['label'],iteration=summary['iteration'],
                        original=summary['true_relative_residual'],canonical=summary['independently_recomputed_true_relative_residual'],
                        fortran_layout=fortran,extended=extended,matrix_sha256=summary['matrix_sha256']))
        result['paths'][folder.name]=dict(all_system_summaries_aligned_to_attempts=total,accepted_attempt_systems=accepted,
            accepted_original_gate_exceedances=original_over,accepted_canonical_gate_exceedances=canonical_over,
            maximum_accepted_original_residual=maximum_original,maximum_accepted_canonical_residual=maximum_canonical,
            retained_actual_systems=len(retained),retained_fortran_evaluation_exact_original_matches=f_matches,
            retained_accepted_actual_systems=retained_accepted,retained_accepted_extended_gate_exceedances=extended_over,
            maximum_retained_accepted_extended_residual=extended_maximum,retained_accepted_canonical_exceedance_examples=examples)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.write_text(json.dumps(audit(a.evidence.resolve()),indent=2,allow_nan=False)+'\n')
