"""Compare actual partial-checkpoint restart fields and complete stored systems."""
import gzip
import hashlib
import json


def rows(path):
    with gzip.open(path,'rt',encoding='utf-8') as stream:
        for line in stream:yield json.loads(line)


def audit(source):
    direct=source/'native-long';resumed=source/'native-long-resumed'
    checkpoint=json.loads((direct/'middle.json').read_text(encoding='utf-8'))
    first=checkpoint['payload']['step']+1
    counts={}
    for name in ['samples','events','linear-systems','tangents']:
        left=(row for row in rows(direct/(name+'.jsonl.gz')) if row['step']>=first)
        right=rows(resumed/(name+'.jsonl.gz'))
        count=0
        while True:
            a=next(left,None);b=next(right,None)
            if a is None or b is None:
                if a is not None or b is not None:raise ValueError('partial restart stream count differs: '+name)
                break
            # Attempt ordinals start at 1 in a new process; numerical contents,
            # accepted steps, Newton iteration numbers and actual systems stay.
            a.pop('attempt',None);b.pop('attempt',None)
            if a!=b:raise ValueError('partial restart stream differs: '+name)
            count+=1
        counts[name]=count
    if not counts['samples'] or not counts['linear-systems']:raise ValueError('empty restart comparison')
    a=(direct/'restart.json').read_bytes();b=(resumed/'restart.json').read_bytes()
    if a!=b:raise ValueError('final sealed partial checkpoint differs')
    return dict(actual_middle_step=first-1,remaining_stream_rows=counts,entire_saved_fields_exact=True,
                actual_linear_systems_exact=True,final_checkpoint_bytes_exact=True,
                final_checkpoint_sha256=hashlib.sha256(a).hexdigest(),
                attempt_ordinals_ignored=True,restart_to_completion_verified=False,numerical_acceptance=False)
