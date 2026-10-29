"""Compare every Richards report event with only its declared timer excluded.

The frozen Rust resume validator normalizes this exact field. No head, flux,
storage, surface decision, attempt, input, physical time or runtime is excluded.
The reference must be a separately hashed and independently audited report.
"""
from pathlib import Path
from itertools import zip_longest
import hashlib,math
import ijson
from streamed_stage_source import CheckedReader

MAX_SOURCE=32*1024**3
TIMER='attempts.item.iterations.item.linear.elapsed_seconds'

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def compare(left,right,*,expected_left_sha256):
    left,right=Path(left),Path(right)
    binding={str(q.resolve()):sha(q) for q in (left,right)}
    if binding[str(left.resolve())]!=expected_left_sha256:raise ValueError('audited reference hash differs')
    if any(q.stat().st_size>MAX_SOURCE for q in (left,right)):raise ValueError('report exceeds source byte budget')
    missing=object();events=0;timers=0;frames=0;numeric_values=0;complete=None;schema=None;last_time=None
    with left.open('rb') as a,right.open('rb') as b:
        streams=[ijson.parse(CheckedReader(f,MAX_SOURCE),use_float=True) for f in (a,b)]
        try:
            for x,y in zip_longest(*streams,fillvalue=missing):
                if x is missing or y is missing:raise ValueError('report event count differs')
                prefix,event,value=x
                if (prefix,event)!=(y[0],y[1]):raise ValueError('report structure differs at '+prefix)
                if event=='number':
                    if not math.isfinite(value) or not math.isfinite(y[2]):raise ValueError('nonfinite report number')
                    numeric_values+=1
                if prefix==TIMER and event=='number':
                    if min(value,y[2])<0:raise ValueError('negative solver timer')
                    timers+=1
                elif value!=y[2]:raise ValueError('non-time report field differs at '+prefix)
                if prefix=='frames.item' and event=='start_map':frames+=1
                if prefix=='complete' and event=='boolean':complete=value
                if prefix=='schema' and event=='string':schema=value
                if prefix=='last_committed.time_s' and event=='number':last_time=value
                events+=1
        except ijson.JSONError as error:raise ValueError('invalid report JSON') from error
    if binding!={str(q.resolve()):sha(q) for q in (left,right)}:raise ValueError('report changed during comparison')
    return dict(passed=True,events=events,numeric_values=numeric_values,frames=frames,
        ignored_solver_timer_values=timers,only_ignored_prefix=TIMER,
        schema=schema,complete=complete,last_committed_time_s=last_time,
        all_non_time_events_equal=True,whole_reports_not_materialized=True,
        source_sha256=binding,source_bytes={str(q.resolve()):q.stat().st_size for q in (left,right)})
