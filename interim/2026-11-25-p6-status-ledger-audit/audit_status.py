#!/usr/bin/env python3
"""Read-only P6 evidence audit. Standard library only; never runs a solver."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

BASE = '79d839e3f89331a82d976efc1f23cbe2c848cbee'
DEST = 'interim/2026-11-25-p6-status-ledger-audit'
STATES = ('限定受入', '未達', '未対応（入力拒否）', '判定不能', '未確認')
DOCS = ['docs/P6_STATUS.ja.md', 'docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md',
        'interim/2026-11-16-issue-triage/todo.ja.md']
REQUIRED = ({f'P3b-{i:02}' for i in range(1, 6)} |
            {'P4a-01', 'P4a-02', 'P4a-03', 'P4b-01', 'P4c-01', 'P4c-02'} |
            {f'STATUS:P6-{i:02}' for i in range(1, 6)} |
            {'NC-01'} | {f'{c}-{i:02}' for c, n in [('A', 4), ('B', 3), ('C', 4)] for i in range(1, n+1)} |
            {f'TODO:P6-{i:02}' for i in range(1, 4)} |
            {'構造・水・JOINT全体', '一般3D空間収束', 'P6全体', 'PC-01', 'CI-01'})


def sha(data):
    return hashlib.sha256(data).hexdigest()


def pointer(obj, key):
    """RFC6901 pointer, including array indexes. Empty pointer denotes root."""
    if key == '':
        return obj
    if not isinstance(key, str) or not key.startswith('/'):
        raise ValueError('not a JSON pointer')
    for token in key[1:].split('/'):
        if re.search(r'~(?![01])', token):
            raise ValueError('invalid pointer escape')
        token = token.replace('~1', '/').replace('~0', '~')
        if isinstance(obj, list):
            if not re.fullmatch(r'0|[1-9][0-9]*', token):
                raise ValueError('invalid array index')
            obj = obj[int(token)]
        else:
            obj = obj[token]
    return obj


def exact(a, b):
    return json.dumps(a, sort_keys=True, ensure_ascii=False) == json.dumps(b, sort_keys=True, ensure_ascii=False)


def safe_path(root, path):
    p = Path(path)
    if p.is_absolute() or '..' in p.parts:
        raise ValueError('unsafe evidence path')
    resolved = (root / p).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError('evidence escapes repository')
    return resolved


def git(root, *args):
    return subprocess.check_output(['git', '-C', str(root), *args], stderr=subprocess.PIPE,
                                   env=dict(os.environ, GIT_NO_LAZY_FETCH='1'))


def audit(root, ledger, *, verify_git=True):
    """Return findings; findings do not change classifications or source bytes."""
    root = Path(root)
    findings = []
    def add(code, message, *, rows=(), evidence=(), severity='error'):
        findings.append(dict(code=code, severity=severity, rows=list(rows), evidence=list(evidence), message=message))
    if ledger.get('base_commit') != BASE:
        add('BASE_COMMIT', '指定起点コミットと異なる')
    inventory = ledger.get('sources', {})
    actual_sources = set(DOCS)
    if (root / 'interim').is_dir():
        actual_sources |= {p.relative_to(root).as_posix() for p in (root/'interim').rglob('*.json')
                           if p.name in ('acceptance.json', 'validation.json', 'measurements.json')
                           and DEST not in p.relative_to(root).as_posix()}
    listed = {s.get('path') for s in inventory.values()}
    for path in sorted(actual_sources - listed):
        add('SOURCE_NOT_INVENTORIED', path)
    cache = {}
    for sid, source in inventory.items():
        path = source.get('path', '')
        try:
            data = safe_path(root, path).read_bytes()
            cache[sid] = data
        except (OSError, ValueError) as e:
            add('SOURCE_MISSING', f'{path}: {e}', evidence=[sid])
            continue
        if sha(data) != source.get('sha256'):
            add('SOURCE_SHA256', path, evidence=[sid])
        if source.get('commit') != BASE:
            add('SOURCE_COMMIT', path, evidence=[sid])
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob != source.get('git_blob_sha1'):
            add('SOURCE_BLOB', path, evidence=[sid])
        if verify_git:
            try:
                fixed = git(root, 'show', f'{BASE}:{path}')
                if fixed != data:
                    add('BASE_BYTES', path, evidence=[sid])
            except subprocess.CalledProcessError:
                add('BASE_EVIDENCE_MISSING', path, evidence=[sid])
        if source.get('kind') == 'json':
            try:
                obj = json.loads(data)
                if isinstance(obj, dict):
                    for key in ('p6_complete', 'P6_complete'):
                        if obj.get(key) is True:
                            add('UNSUPPORTED_COMPLETION', f'{path}#/{key}: P6全体の完了証拠なし', evidence=[sid])
                    if obj.get('passed') is True and any(obj.get(k) is False for k in
                            ('requested_all_combinations_passed', 'requested_all_items_passed')):
                        add('LIMITED_PASS_NOT_OVERALL', f'{path}: passed=true は要求全体の成功ではない',
                            evidence=[sid], severity='info')
                for proof in source.get('status_keys', []):
                    if not exact(pointer(obj, proof['pointer']), proof['value']):
                        add('INVENTORY_KEY_VALUE', f"{path}#{proof['pointer']}", evidence=[sid])
            except (ValueError, KeyError, IndexError, TypeError) as e:
                add('INVENTORY_JSON', f'{path}: {e}', evidence=[sid])
    for sid, data in cache.items():
        if inventory[sid].get('kind') != 'markdown':
            continue
        for lineno, line in enumerate(data.decode('utf-8').splitlines(), 1):
            # Only the status column of an identified acceptance row; negated/historical prose is not a completion claim.
            cells = [c.strip() for c in line.split('|')]
            if len(cells) >= 5 and re.fullmatch(r'P(?:3b|4[abc]|6)-[0-9]+', cells[1]):
                status = cells[3].replace('**', '')
                if re.search(r'(?<!未)(?<!不)完了|全体受入済み|全面合格', status):
                    add('UNSUPPORTED_COMPLETION', f"{inventory[sid]['path']}:{lineno}: {status}", evidence=[sid])
    proofs = ledger.get('evidence', {})
    valid = set()
    for eid, proof in proofs.items():
        sid = proof.get('source')
        if sid not in cache:
            add('EVIDENCE_SOURCE', '根拠ソースが存在しない', evidence=[eid])
            continue
        try:
            data = cache[sid]
            if 'pointer' in proof:
                got = pointer(json.loads(data), proof['pointer'])
            else:
                lines = data.decode('utf-8').splitlines()
                start, end = proof['lines']
                if not (1 <= start <= end <= len(lines)):
                    raise ValueError('line range')
                got = '\n'.join(lines[start-1:end])
            if not exact(got, proof['value']):
                add('EVIDENCE_VALUE', 'JSONキー値または引用行が不一致', evidence=[eid])
            else:
                valid.add(eid)
        except (ValueError, KeyError, IndexError, TypeError) as e:
            add('EVIDENCE_SELECTOR', str(e), evidence=[eid])
    rows = ledger.get('rows', [])
    ids = [r.get('id') for r in rows]
    if len(set(ids)) != len(ids):
        add('DUPLICATE_ROW', '台帳IDの重複')
    targets = {r.get('target') for r in rows}
    for missing in sorted(REQUIRED-targets):
        add('COVERAGE_MISSING', missing)
    grouped = defaultdict(list)
    for row in rows:
        rid, state = row.get('id'), row.get('state')
        if state not in STATES:
            add('INVALID_STATE', str(state), rows=[rid])
        eids = row.get('evidence', [])
        if not eids or any(e not in valid for e in eids):
            add('ROW_WITHOUT_EVIDENCE', '検証可能な行根拠が不足', rows=[rid], evidence=eids)
        if not row.get('scope') or not row.get('reason'):
            add('ROW_SCOPE_MISSING', '範囲または判定理由が欠落', rows=[rid])
        if state == '限定受入':
            supporting = [proofs[e] for e in eids if e in valid]
            accepted = any(p.get('supports') == '限定受入' for p in supporting)
            if not accepted:
                add('UNSUPPORTED_ACCEPTANCE', '限定受入の明示証拠がない', rows=[rid])
            allowed = {p.get('scope') for p in supporting if p.get('supports') == '限定受入'}
            if row.get('scope') not in allowed:
                add('SCOPE_EXPANSION', '限定受入の証拠と行の適用範囲が異なる', rows=[rid])
        if row.get('completed') is not False:
            add('UNSUPPORTED_COMPLETION', 'この台帳で全体完了は認定しない', rows=[rid])
        grouped[(row.get('target'), row.get('scope'))].append(row)
    for (target, scope), group in grouped.items():
        states = {r.get('state') for r in group}
        if '限定受入' in states and states & {'未達', '未対応（入力拒否）'}:
            add('ROW_STATE_CONTRADICTION', f'{target}: {scope}', rows=[r['id'] for r in group])
    # Compare explicit, reviewable document claims. Same-scope disagreement is never silently reconciled.
    claims = ledger.get('claims', [])
    claim_groups = defaultdict(list)
    for claim in claims:
        eid = claim.get('evidence')
        if eid not in valid:
            add('CLAIM_EVIDENCE', '文書主張の引用が不正', evidence=[eid])
            continue
        claim_groups[(claim.get('target'), claim.get('scope'))].append(claim)
        value = proofs[eid]['value']
        if claim.get('state') == '完了' or (isinstance(value, str) and
                re.search(r'(?:^|[。|])(?:完了|全体受入済み)(?:[。|]|$)', value)):
            add('UNSUPPORTED_COMPLETION', '全体完了を支える受入証拠がない', evidence=[eid])
    for (target, scope), group in sorted(claim_groups.items()):
        for i, a in enumerate(group):
            for b in group[i+1:]:
                incompatible = ({a.get('state'), b.get('state')} & {'限定受入', '完了'} and
                                {a.get('state'), b.get('state')} & {'未達', '未対応（入力拒否）'})
                if incompatible:
                    historical = a.get('historical') or b.get('historical')
                    add('HISTORICAL_SCOPE_CONFLICT' if historical else 'DOCUMENT_CONTRADICTION',
                        f'{target}: {scope}。時点・範囲を併読して手動確認。自動修正なし。',
                        evidence=[a['evidence'], b['evidence']], severity='warning' if historical else 'error')
    for observation in ledger.get('document_observations', []):
        eids = observation.get('evidence', [])
        if not eids or any(e not in valid for e in eids):
            add('OBSERVATION_EVIDENCE', '文書間差分の引用が不正', evidence=eids)
        else:
            add(observation['code'], observation['message'], evidence=eids, severity='warning')
    external = ledger.get('external_references', {})
    for key in ('ci_snapshot', 'parallel_prs'):
        if key in external:
            try:
                data = safe_path(root, DEST+'/'+external[key]).read_bytes()
                if sha(data) != external[key+'_sha256']:
                    add('EXTERNAL_SNAPSHOT_SHA', external[key])
                if key == 'ci_snapshot' and not exact(json.loads(data)['workflow_runs'], ledger.get('ci_records', [])):
                    add('CI_SNAPSHOT_VALUE', 'CI記録が保存スナップショットと不一致')
            except (OSError, ValueError, KeyError) as e:
                add('EXTERNAL_SNAPSHOT', str(e))
    # Any positive CI claim requires an exact successful, completed run at the claimed head.
    for ci in ledger.get('ci_records', []):
        if ci.get('claimed_success'):
            if ci.get('status') != 'completed' or ci.get('conclusion') != 'success':
                add('CI_NOT_SUCCESS', f"run {ci.get('id')}: {ci.get('status')}/{ci.get('conclusion')}")
            if not re.fullmatch(r'[a-f0-9]{40}', str(ci.get('head_sha', ''))) or ci.get('head_sha') != ci.get('claimed_head_sha'):
                add('CI_WRONG_HEAD', f"run {ci.get('id')}")
        if ci.get('conclusion') in ('cancelled', 'timed_out', 'failure'):
            add('CI_NON_SUCCESS_RECORD', f"run {ci.get('id')}: {ci['conclusion']}。成功に数えない。", severity='info')
    if verify_git:
        try:
            diff = git(root, 'diff', '--cached', '--no-renames', '--name-status', BASE, '--').decode()
            modified = git(root, 'ls-files', '--modified').decode()
            for path in modified.splitlines():
                if not path.startswith(DEST+'/'):
                    add('ADDITIONS_ONLY', 'modified tracked file: '+path)
            untracked = git(root, 'ls-files', '--others', '--exclude-standard').decode()
            for path in untracked.splitlines():
                if not path.startswith(DEST+'/'):
                    add('ADDITIONS_ONLY', 'untracked outside destination: '+path)
            for line in diff.splitlines():
                kind, path = line.split('\t', 1)
                if kind != 'A' or not path.startswith(DEST+'/'):
                    add('ADDITIONS_ONLY', line)
        except (subprocess.CalledProcessError, ValueError) as e:
            add('DIFF_UNAVAILABLE', str(e))
    counts = Counter(r.get('state') for r in rows)
    return {'base_commit': BASE, 'row_count': len(rows),
            'state_counts': {s: counts[s] for s in STATES},
            'source_count': len(inventory), 'evidence_count': len(proofs),
            'finding_counts': dict(Counter(f['code'] for f in findings)),
            'integrity_ok': not any(f['severity'] == 'error' for f in findings),
            'numerical_acceptance_performed': False, 'findings': findings}


def table(ledger):
    def cell(value):
        return str(value).replace('|', '\\|').replace('\n', '<br>')
    lines = ['# P6 状態台帳（日本語表）', '',
             '起点: `'+ledger['base_commit']+'`。限定受入は既存の限定範囲のみ。全体完了の認定はしない。', '',
             '|行ID / 対象|状態|範囲|判定理由・残り|根拠ID|', '|---|---|---|---|---|']
    for r in ledger['rows']:
        lines.append('|'+ '|'.join(cell(v) for v in [r['id']+' / '+r['target'], r['state'], r['scope'], r['reason'], ', '.join(r['evidence'])])+'|')
    lines += ['', '根拠IDは ledger.json の evidence を参照。source にパス・起点コミット・SHA-256・Git blob SHA-1、pointer にJSONキー、lines に原文行範囲を保持する。', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--ledger', type=Path, default=Path(__file__).with_name('ledger.json'))
    parser.add_argument('--table', action='store_true', help='日本語表を標準出力へ出す（書込なし）')
    args = parser.parse_args()
    try:
        ledger = json.loads(args.ledger.read_text(encoding='utf-8'))
        if args.table:
            print(table(ledger), end='')
            return 0
        result = audit(args.root, ledger)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0 if result['integrity_ok'] else 1
    except (OSError, ValueError, TypeError, KeyError) as e:
        print(json.dumps({'integrity_ok': False, 'error': str(e)}, ensure_ascii=False))
        return 2


if __name__ == '__main__':
    sys.exit(main())
