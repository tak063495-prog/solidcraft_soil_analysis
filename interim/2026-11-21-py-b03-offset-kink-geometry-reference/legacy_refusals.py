"""Read archived B-03 inputs; independently check geometry, without solver calls."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

from geometry_reference import dump, orientation_diagnostics

BASE = 'interim/2026-11-14-plate6-curved-offset'
PREFIX = 'retained-failures/scope-refusals/'
CASES = ('outside_half_thickness', 'folded_offset_normals')


def inspect_direct_geometry(stages):
    cfg = stages['beam_system']
    positions = {r['id']: np.asarray(r['position_m'], float) for r in cfg['nodes']}
    sections = {r['id']: r for r in cfg['plate_sections']}
    faces, nodes, errors = [], {}, []
    for plate in cfg['plates']:
        if plate['formulation'] != 'mindlin_tri6_offset':
            errors.append('outside_explicit_offset_formulation')
            continue
        ids = plate['nodes']
        if len(ids) != 6:
            raise ValueError('only archived direct six-node declarations are interpreted')
        p = np.array([positions[i] for i in ids])
        n = np.cross(p[1]-p[0], p[2]-p[0])
        if not np.isfinite(p).all() or not np.isfinite(n).all() or np.linalg.norm(n) == 0:
            raise ValueError('invalid/degenerate archived face')
        n /= np.linalg.norm(n)
        e = plate.get('attachment_offset_m', 0.)
        t = sections[plate['section_id']]['thickness_m']
        face_errors = orientation_diagnostics(n, e, t, n)
        errors.extend(face_errors)
        en = e*n
        faces.append({'face_id': plate['id'], 'normal': n.tolist(), 'e_m': e, 't_m': t,
                      'signed_normal_offset_m': en.tolist(), 'errors': face_errors})
        for node, point in zip(ids, p):
            nodes.setdefault(node, []).append((plate['id'], en, point-en, n))
    shared = []
    for node, uses in nodes.items():
        if len(uses) < 2:
            continue
        first = uses[0]
        for other in uses[1:]:
            delta = other[1]-first[1]
            inconsistent = bool(np.any(delta != 0))
            if inconsistent:
                errors.append('shared_signed_normal_offset_mismatch')
            shared.append({'shared_node_id': node, 'face_ids': [first[0], other[0]],
                           'signed_offset_delta_m': delta.tolist(),
                           'neutral_point_delta_m': (other[2]-first[2]).tolist(),
                           'normal_dot': float(first[3]@other[3]),
                           'inconsistent': inconsistent})
    return {'rejected': bool(errors), 'errors': sorted(set(errors)), 'faces': faces,
            'shared_nodes': shared, 'accepted': None,
            'scope': 'independent exact geometry guard; not a complete production input validator'}


def audit_repository(root):
    root = Path(root)
    source = root/BASE
    manifest = json.loads((source/'evidence-manifest.json').read_text(encoding='utf-8'))
    expected = {r['path']: r for r in manifest['files']}
    old = json.loads((source/'retained-failures/scope-rejections.json').read_text(encoding='utf-8'))
    records = {r['case']: r for r in old['records']}
    archive = source/'evidence.zip'
    reports = []
    with zipfile.ZipFile(archive) as z:
        for name in CASES:
            checked = []
            documents = {}
            for field in ('mesh', 'ground', 'stages'):
                member = PREFIX+name+'/'+field+'.json'
                raw = z.read(member)
                digest = hashlib.sha256(raw).hexdigest()
                if digest != expected[member]['sha256'] or len(raw) != expected[member]['bytes']:
                    raise ValueError('archived input hash/length mismatch: '+member)
                documents[field] = json.loads(raw)
                checked.append({'archive_member': member, 'sha256': digest, 'bytes': len(raw)})
            report = inspect_direct_geometry(documents['stages'])
            record = records[name]
            reports.append({'case': name, 'input_status': '確認済み', 'inputs': checked,
                            'independent_guard': report,
                            'historical_body_refusal': record,
                            'same_rejection_confirmed': bool(report['rejected'] and record['refused'] and record['returncode'] != 0),
                            'body_rerun_by_this_script': False})
    return {'accepted': None, 'b03_complete': False,
            'source_base': BASE, 'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
            'manifest_sha256': hashlib.sha256((source/'evidence-manifest.json').read_bytes()).hexdigest(),
            'cases': reports,
            'unconfirmed': [
                {'item': 'noncoplanar physical kink input', 'status': '未確認',
                 'reason': 'folded_offset_normals is a coplanar winding reversal, not a finite-angle kink'},
                {'item': 'parallel physical step input', 'status': '未確認',
                 'reason': 'no separately identified B-03 step input in this cited archive'},
                {'item': 'compact/Native/GUI/GPU/real CAD input and display', 'status': '未確認',
                 'reason': 'PC procedure only; no fabricated historical input'}]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository_root')
    parser.add_argument('output')
    args = parser.parse_args()
    try:
        report = audit_repository(args.repository_root)
        dump(report, args.output)
        return 2 if all(r['same_rejection_confirmed'] for r in report['cases']) else 1
    except (KeyError, ValueError, OSError, zipfile.BadZipFile) as exc:
        dump({'accepted': None, 'status': '未確認', 'error': str(exc)}, args.output)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
