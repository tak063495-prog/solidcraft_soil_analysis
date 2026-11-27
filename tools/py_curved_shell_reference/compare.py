"""Compare a future solver's nodal CSV/JSON against a fixed literature probe.

Exit 0: passes a documented, defined tolerance; 1: fails it; 2: undefined;
3: invalid/missing/ambiguous input. No user-supplied tolerance override.
CSV columns: node_id(optional),x,y,z,ux,uy,uz. JSON: a list of equivalent
objects or {"benchmark_id": "...", "nodes": [...]}. Coordinates and vector
components use the global frame and consistent units in benchmarks.json.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

from reference import catalog


COLUMNS = ('x','y','z','ux','uy','uz')


def read_nodes(path, benchmark_id):
    path=Path(path)
    if path.suffix.lower()=='.csv':
        with path.open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.DictReader(f))
    elif path.suffix.lower()=='.json':
        obj=json.loads(path.read_text(encoding='utf-8-sig'))
        if isinstance(obj,dict):
            if obj.get('benchmark_id',benchmark_id)!=benchmark_id:
                raise ValueError('benchmark_id mismatch')
            rows=obj['nodes']
        else:
            rows=obj
    else:
        raise ValueError('expected .csv or .json')
    if not isinstance(rows,list) or not rows:
        raise ValueError('nodes must be a nonempty list')
    data=[];ids=set();points=set()
    for row in rows:
        if not isinstance(row,dict) or any(c not in row for c in COLUMNS):
            raise ValueError('each node requires x,y,z,ux,uy,uz')
        if any(isinstance(row[c],bool) for c in COLUMNS):
            raise ValueError('boolean is not a nodal number')
        values=[float(row[c]) for c in COLUMNS]
        if not np.isfinite(values).all():
            raise ValueError('nonfinite nodal value')
        point=tuple(values[:3])
        if point in points:
            raise ValueError('duplicate coordinates')
        points.add(point)
        if 'node_id' in row:
            node_id=str(row['node_id'])
            if not node_id or node_id in ids:
                raise ValueError('empty or duplicate node_id')
            ids.add(node_id)
        data.append(values)
    return np.asarray(data)


def compare(case, data, acceptance):
    probe=case['probe']
    matching=np.flatnonzero(np.linalg.norm(data[:,:3]-probe['xyz'],axis=1)<=probe['coordinate_tolerance'])
    if len(matching)!=1:
        raise ValueError('probe must match exactly one physical node; no nearest-node extrapolation')
    value=float(data[matching[0],COLUMNS.index(probe['component'])])
    target=case['reference']['signed_value']
    absolute=abs(value-target)
    relative=absolute/abs(target)
    tolerance=acceptance['relative_displacement_tolerance']
    passed=None if tolerance is None else bool(relative<=tolerance)
    return dict(benchmark_id=case['id'],probe=probe,computed=value,reference=target,
                signed_difference=value-target,absolute_difference=absolute,
                relative_difference=relative,relative_tolerance=tolerance,
                status='未定義' if passed is None else ('合格' if passed else '不合格'),
                passed=passed,B_01_passed=False,
                tolerance_basis=acceptance,
                reference_sources=case['reference']['source_ids'],
                scope='single literature displacement probe only; not a full-field or B-01 acceptance')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--benchmark',required=True,choices=[c['id'] for c in catalog()['benchmarks']])
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--out',type=Path)
    args=p.parse_args(argv)
    try:
        meta=catalog()
        case=next(c for c in meta['benchmarks'] if c['id']==args.benchmark)
        result=compare(case,read_nodes(args.input,args.benchmark),meta['acceptance'])
        code=2 if result['passed'] is None else (0 if result['passed'] else 1)
    except (ValueError,TypeError,KeyError,OSError) as e:
        result=dict(status='入力不正',passed=False,B_01_passed=False,error=str(e));code=3
    text=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.write_text(text,encoding='utf-8')
    else:
        print(text,end='')
    return code


if __name__=='__main__':
    sys.exit(main())
