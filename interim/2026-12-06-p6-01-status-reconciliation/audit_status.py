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
import tomllib

BASE = '8342762bd714e33523ce5fef772572699fc21881'
DEST = 'interim/2026-12-06-p6-01-status-reconciliation'
SUPPLEMENTAL_PATH = 'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'
SUPPLEMENTAL_SHA256 = '1824329ba9422412c7314e657de014ac21bfbf04952c96a41a3b62823d134184'
PHYSICAL_APPENDICES_PATH = 'interim/2026-11-13-nc-lambda-physical-limit/document-appendices.json'
PHYSICAL_APPENDICES_SHA256 = '2984747e89a2d3e2ff3f45ae3b32c06c0d0522a85e844b4ca349705ea1736bdc'
COLLAPSE_APPENDICES_PATH = 'interim/2026-11-13-nc-lambda-collapse-classification/document-appendices.json'
COLLAPSE_APPENDICES_SHA256 = '9ef6ecf6de7eafb6dc199172249f62179812736c8c03b7c780a1dc1eed3f2bfc'
STATES = ('限定受入', '未達', '未対応（入力拒否）', '判定不能', '未確認・PC手順待ち')
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


def collapse_document_appendices(root):
    """Only the exact pinned new suffix can precede historical hash checks."""
    root = Path(root)
    path = root / COLLAPSE_APPENDICES_PATH
    if not path.parent.exists():
        return {}
    raw = path.read_bytes()
    if sha(raw) != COLLAPSE_APPENDICES_SHA256:
        raise ValueError('NC collapse document appendix manifest SHA-256 differs')
    data = json.loads(raw)
    if (set(data) != {'schema', 'diagnostic_only', 'passed', 'nc01_complete', 'appendices'}
            or data['schema'] != 'nc-collapse-document-appendices/1'
            or data['diagnostic_only'] is not True or data['passed'] is not False
            or data['nc01_complete'] is not False):
        raise ValueError('NC collapse appendix status or schema differs')
    tails = data['appendices']
    if (set(tails) != {'docs/LOW_PRELOAD_NC.ja.md', 'interim/2026-11-16-issue-triage/todo.ja.md'}
            or any(not isinstance(t, str) or not t for t in tails.values())):
        raise ValueError('NC collapse appendix scope differs')
    for name, tail in tails.items():
        if not safe_path(root, name).read_bytes().endswith(tail.encode('utf-8')):
            raise ValueError('NC collapse document appendix suffix differs: '+name)
    return tails


def physical_document_appendices(root):
    """Validate a new, pinned tail without resealing any historical document."""
    root = Path(root)
    raw = (root / PHYSICAL_APPENDICES_PATH).read_bytes()
    if sha(raw) != PHYSICAL_APPENDICES_SHA256:
        raise ValueError('NC physical document appendix manifest SHA-256 differs')
    data = json.loads(raw)
    if (set(data) != {'schema', 'diagnostic_only', 'passed', 'nc01_complete', 'appendices'}
            or data['schema'] != 'nc-physical-document-appendices/1'
            or data['diagnostic_only'] is not True or data['passed'] is not False
            or data['nc01_complete'] is not False):
        raise ValueError('NC physical appendix status or schema differs')
    tails = data['appendices']
    if (set(tails) != {'docs/LOW_PRELOAD_NC.ja.md', 'interim/2026-11-16-issue-triage/todo.ja.md'}
            or any(not isinstance(t, str) or not t for t in tails.values())):
        raise ValueError('NC physical appendix scope differs')
    collapse = collapse_document_appendices(root)
    combined = {name: tail+collapse.get(name, '') for name, tail in tails.items()}
    for name, tail in combined.items():
        if not safe_path(root, name).read_bytes().endswith(tail.encode('utf-8')):
            raise ValueError('NC physical document appendix suffix differs: '+name)
    return combined


def audit(root, ledger, *, verify_git=True):
    """Return findings; findings do not change classifications or source bytes."""
    root = Path(root)
    findings = []
    def add(code, message, *, rows=(), evidence=(), severity='error'):
        findings.append(dict(code=code, severity=severity, rows=list(rows), evidence=list(evidence), message=message))
    if ledger.get('base_commit') != BASE:
        add('BASE_COMMIT', '指定起点コミットと異なる')
    inventory = dict(ledger.get('sources', {}))
    # Append new snapshots without resealing the historical BASE ledger. The
    # manifest is pinned here; each added JSON still uses every source check.
    supplemental_ids = set()
    manifest = root / SUPPLEMENTAL_PATH
    if manifest.exists() or (root / DEST / 'ledger.json').exists():
        try:
            raw = manifest.read_bytes()
            if sha(raw) != SUPPLEMENTAL_SHA256:
                raise ValueError('supplemental inventory SHA-256 differs')
            extension = tomllib.loads(raw.decode('utf-8'))
            if extension.get('schema') != 'p601-inventory/v1':
                raise ValueError('unsupported supplemental inventory schema')
            paths = {s.get('path') for s in inventory.values()}
            pending = {}
            for source in extension['sources']:
                source = dict(source)
                sid = source.pop('id')
                path = source.get('path', '')
                if (set(source) != {'path', 'commit', 'kind', 'sha256', 'git_blob_sha1', 'status_keys'}
                        or not isinstance(sid, str) or not sid or sid in inventory or sid in pending
                        or path in paths or not path.startswith('interim/')
                        or source.get('kind') != 'json' or source.get('origin') is not None
                        or not re.fullmatch(r'[a-f0-9]{40}', source.get('commit', ''))
                        or not re.fullmatch(r'[a-f0-9]{64}', source.get('sha256', ''))
                        or not re.fullmatch(r'[a-f0-9]{40}', source.get('git_blob_sha1', ''))
                        or not source.get('status_keys')):
                    raise ValueError('invalid or duplicate supplemental source: '+str(sid))
                safe_path(root, path)
                for proof in source['status_keys']:
                    proof['value'] = json.loads(proof.pop('value_json'))
                pending[sid] = source
                paths.add(path)
            if not pending:
                raise ValueError('empty supplemental inventory')
            inventory.update(pending)
            supplemental_ids = set(pending)
        except (OSError, ValueError, KeyError, TypeError) as e:
            add('SUPPLEMENTAL_INVENTORY', str(e))
    actual_sources = set(DOCS)
    if (root / 'interim').is_dir():
        actual_sources |= {p.relative_to(root).as_posix() for p in (root/'interim').rglob('*.json')
                           if p.name in ('acceptance.json', 'validation.json', 'measurements.json')
                           and DEST not in p.relative_to(root).as_posix()}
    listed = {s.get('path') for s in inventory.values()}
    for path in sorted(actual_sources - listed):
        add('SOURCE_NOT_INVENTORIED', path)
    # Later append-only reconciliation; the historical ledger stays byte-identical.
    later = {}
    tail_manifest = root / 'interim/2026-12-07-p6-03-cloud-reconciliation/document-appendices.json'
    if tail_manifest.exists():
        raw_tail = tail_manifest.read_bytes()
        if sha(raw_tail) != 'e82a012789159ef513706ef1377f31f0e5c3548a90aed4ce02f605c3f50207fe':
            add('LATER_APPENDIX_MANIFEST', 'P6-03追補manifestの厳密SHAが異なる')
        else:
            later = json.loads(raw_tail)
            if set(later) != {DOCS[0], DOCS[2]} or any(not isinstance(t, str) or not t for t in later.values()):
                add('LATER_APPENDIX_SCOPE', 'P6-03追補の対象/内容が異なる')
                later = {}
    physical_tails = {}
    if (root / PHYSICAL_APPENDICES_PATH).exists() or (root / DEST / 'ledger.json').exists():
        try:
            physical_tails = physical_document_appendices(root)
        except (OSError, ValueError, KeyError, TypeError) as e:
            add('PHYSICAL_APPENDIX_MANIFEST', str(e))
    cache = {}
    for sid, source in inventory.items():
        path = source.get('path', '')
        if source.get('origin') == 'github_snapshot' and (path != DEST+'/ci-snapshot.json' or source.get('commit') != 'e278dbf17c771e87098b6180e9c65d9e4a6ad2e7'):
            add('EXTERNAL_SOURCE_SCOPE', '固定Git証拠を外部snapshotへ再分類できない', evidence=[sid])
        try:
            data = safe_path(root, path).read_bytes()
            if path in physical_tails:
                # The complete tail was checked above, before old-byte checks.
                data = data[:-len(physical_tails[path].encode('utf-8'))]
            if path in later:
                tail = later[path].encode('utf-8')
                if not data.endswith(tail):
                    add('LATER_DOCUMENT_APPENDIX', path, evidence=[sid])
                else:
                    data = data[:-len(tail)]
            if source.get('appendix') is not None:
                appendix = source['appendix'].encode('utf-8')
                if not data.endswith(appendix):
                    add('DOCUMENT_APPENDIX', path, evidence=[sid])
                else:
                    data = data[:-len(appendix)]
            cache[sid] = data
        except (OSError, ValueError) as e:
            add('SOURCE_MISSING', f'{path}: {e}', evidence=[sid])
            continue
        if sha(data) != source.get('sha256'):
            add('SOURCE_SHA256', path, evidence=[sid])
        fixed_commit = source['commit'] if sid in supplemental_ids else BASE
        if source.get('origin') != 'github_snapshot' and source.get('commit') != fixed_commit:
            add('SOURCE_COMMIT', path, evidence=[sid])
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if blob != source.get('git_blob_sha1'):
            add('SOURCE_BLOB', path, evidence=[sid])
        if verify_git and source.get('origin') != 'github_snapshot':
            try:
                fixed = git(root, 'rev-parse', '--verify', f'{fixed_commit}:{path}').decode().strip()
                if fixed != blob:
                    add('BASE_BYTES', path, evidence=[sid])
            except subprocess.CalledProcessError:
                add('BASE_EVIDENCE_MISSING', path, evidence=[sid])
        if source.get('kind') == 'json':
            try:
                obj = json.loads(data)
                if isinstance(obj, dict):
                    for key in ('p6_complete', 'P6_complete', 'nc01_complete', 'a01_complete', 'a02_complete', 'a03_complete', 'a04_complete', 'b01_complete', 'b02_complete', 'b03_complete', 'c01_complete', 'c02_complete', 'c03_complete', 'c04_complete'):
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
        complete_document = (data.decode('utf-8') + inventory[sid].get('appendix', '')
                             + later.get(inventory[sid]['path'], '')
                             + physical_tails.get(inventory[sid]['path'], ''))
        for lineno, line in enumerate(complete_document.splitlines(), 1):
            # Only the status column of an identified acceptance row; negated/historical prose is not a completion claim.
            if re.search(r'(?:NC-01|P6全体|P3b/P4各行)\s*(?:は|が|[:：])\s*(?:全て|すべて)?(?:完了(?:しました|です|済み)|全体受入済み|全面合格)', line) and not re.search(r'しない|ではない|認定せず|禁止|引用|誤記', line):
                add('UNSUPPORTED_COMPLETION', f"{inventory[sid]['path']}:{lineno}: 全体完了の肯定表記", evidence=[sid])
            cells = [c.strip() for c in line.split('|')]
            if len(cells) >= 5 and re.fullmatch(r'P(?:3b|4[abc]|6)-[0-9]+', cells[1]):
                status = cells[3].replace('**', '')
                if re.match(r'(?:完了|全体受入済み|全面合格)(?:$|[（(、。・ ])', status):
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
        if row.get('target') in {'NC-01', 'P6全体', '一般3D空間収束', 'TODO:P6-01'} and row.get('scope_kind') == 'overall' and state == '限定受入':
            add('SCOPE_EXPANSION', '診断または限定証拠を全体要件の受入に転用', rows=[rid])
        if state == '限定受入':
            supporting = [proofs[e] for e in eids if e in valid]
            accepted = any(p.get('supports') == '限定受入' and p.get('value') is not False and p.get('value') is not None for p in supporting)
            for proof in supporting:
                if proof.get('supports') == '限定受入' and proof.get('value') in (False, None):
                    add('FALSE_OR_NULL_ACCEPTANCE', 'false/null根拠を限定受入へ転用', rows=[rid])
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
    if 'ci_snapshot' in external:
        critical = [c for c in ledger.get('ci_records', []) if c.get('id') == 43191750591]
        if len(critical) != 1 or critical[0].get('conclusion') != 'failure' or not critical[0].get('infrastructure_failure'):
            add('KNOWN_CI_FAILURE_LOST', 'e278dbf / run43191750591のcheckout失敗事実が失われた')
    # Any positive CI claim requires an exact successful, completed run at the claimed head.
    for ci in ledger.get('ci_records', []):
        if ci.get('claimed_success'):
            jobs = ci.get('jobs', [])
            if ci.get('infrastructure_failure') or any(j.get('conclusion') != 'success' for j in jobs):
                add('CI_INFRA_OR_JOB_FAILURE', f"run {ci.get('id')}: 全job成功を満たさない")
            if ci.get('status') != 'completed' or ci.get('conclusion') != 'success':
                add('CI_NOT_SUCCESS', f"run {ci.get('id')}: {ci.get('status')}/{ci.get('conclusion')}")
            if not re.fullmatch(r'[a-f0-9]{40}', str(ci.get('head_sha', ''))) or ci.get('head_sha') != ci.get('claimed_head_sha'):
                add('CI_WRONG_HEAD', f"run {ci.get('id')}")
        if ci.get('conclusion') in ('cancelled', 'timed_out', 'failure', 'action_required', 'stale'):
            add('CI_NON_SUCCESS_RECORD', f"run {ci.get('id')}: {ci['conclusion']}。成功に数えない。", severity='info')
    counts = Counter(r.get('state') for r in rows)
    return {'base_commit': BASE, 'row_count': len(rows),
            'state_counts': {s: counts[s] for s in STATES},
            'source_count': len(inventory), 'evidence_count': len(proofs),
            'finding_counts': dict(Counter(f['code'] for f in findings)),
            'integrity_ok': not any(f['severity'] == 'error' for f in findings),
            'numerical_acceptance_performed': False, 'findings': findings}


def citation(ledger, eid):
    p = ledger['evidence'][eid]
    s = ledger['sources'][p['source']]
    selector = p.get('pointer') if 'pointer' in p else ':'.join(map(str, p['lines']))
    return f"{eid}: `{s['path']}` #{selector}; commit `{s['commit']}`; SHA256 `{s['sha256']}`; blob `{s['git_blob_sha1']}`"


def table(ledger):
    def cell(value):
        return str(value).replace('|', '\\|').replace('\n', '<br>')
    lines = ['# P6 状態台帳（日本語表）', '',
             '起点: `'+ledger['base_commit']+'`。限定受入は既存の限定範囲のみ。全体完了の認定はしない。', '',
             '|行ID / 対象|状態|範囲|判定理由・残り|根拠ID|', '|---|---|---|---|---|']
    for r in ledger['rows']:
        lines.append('|'+ '|'.join(cell(v) for v in [r['id']+' / '+r['target'], r['state'], r['scope'], r['reason'], '<br>'.join(citation(ledger, e) for e in r['evidence'])])+'|')
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
