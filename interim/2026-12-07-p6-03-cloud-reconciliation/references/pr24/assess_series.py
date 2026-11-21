"""Guarded Richardson diagnosis of future common-sample spatial/load/material series."""
import argparse
import json
from pathlib import Path
from compare_history import read_json
from mcc_reference import richardson


def assess(data):
    axis = data['axis']
    if axis not in ('spatial', 'load_increment', 'material_increment'):
        raise ValueError('unknown refinement axis')
    parameter_unit=data['parameter_unit']
    if parameter_unit not in (('m',) if axis=='spatial' else ('1','kPa')):
        raise ValueError('explicit supported refinement parameter unit required')
    if data['run_status'] != 'completed':
        raise ValueError('incomplete/timeout/cancelled series')
    levels = data['levels']
    if len(levels) != 3:
        raise ValueError('exactly three refinement levels required')
    required = ('physics_sha256', 'observable', 'unit', 'sample_ids', 'depth_m',
                'boundary_id', 'load_path_id', 'load_increment', 'material_increment')
    varying = {'spatial': 'h_m', 'load_increment': 'load_increment',
               'material_increment': 'material_increment'}[axis]
    for k in required:
        if k == varying:
            continue
        v = levels[0][k]
        if v is None or any(level[k] != v for level in levels[1:]):
            raise ValueError('confounded/missing fixed condition: '+k)
    if axis != 'spatial':
        if any(level['mesh_sha256'] != levels[0]['mesh_sha256'] for level in levels):
            raise ValueError('mesh changed in increment series')
    else:
        if data.get('nested_h_refinement') is not True:
            raise ValueError('spatial series must explicitly declare nested in-plane h refinement')
    ids = levels[0]['sample_ids']
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('invalid common sample IDs')
    for level in levels:
        if len(level['values']) != len(ids):
            raise ValueError('missing common samples; no interpolation')
    h = [float(level[varying]) for level in levels]
    if not h[0] > h[1] > h[2] > 0:
        raise ValueError('refinement must decrease')
    r = h[0]/h[1]
    if abs(r-h[1]/h[2]) > 1e-12*r:
        raise ValueError('nonuniform refinement ratio')
    result = richardson([level['values'] for level in levels], ratio=r)
    result.update(axis=axis, parameter_unit=parameter_unit, observable=levels[0]['observable'], unit=levels[0]['unit'],
                  status='diagnostic_only', acceptance_scope='未定義; no new acceptance')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    try:
        result = assess(read_json(a.input))
        with a.output.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write('\n')
    except (ValueError, KeyError, TypeError, OverflowError) as exc:
        p.exit(2, 'REJECTED: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
