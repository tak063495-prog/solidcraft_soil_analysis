"""Read-only A-01 archive/source scope inventory. Never execute archived solvers."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from reference import dump

PINNED = '7ae0a2a669e0e1df93f3dd9791bbeda056334c3e'
FOLDER = Path('interim/2026-11-13-joint-ext-combos')


def audit(root, archive):
    root = Path(root); archive = Path(archive)
    acceptance = json.loads((root/FOLDER/'acceptance.json').read_text(encoding='utf-8'))
    manifest = json.loads((root/FOLDER/'evidence-manifest.json').read_text(encoding='utf-8'))
    artifact = json.loads((root/FOLDER/'artifact-manifest.json').read_text(encoding='utf-8'))
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != artifact['sha256']['evidence.zip']:
        raise ValueError('old archive SHA differs; do not regenerate/replace old evidence')
    entries = {}; logs = []; targeted = []
    with zipfile.ZipFile(archive) as z:
        for name, expected in manifest['sha256']['evidence.zip'].items():
            if not name.startswith('scope-rejections/'):
                continue
            raw = z.read(name)
            got = hashlib.sha256(raw).hexdigest()
            if got != expected:
                raise ValueError('old archive entry differs: '+name)
            entries[name] = got
            if name.endswith('/refused.log'):
                logs.append(dict(path=name, text=raw.decode('utf-8'), sha256=got))
            elif name.endswith('/case.json'):
                case = json.loads(raw)
                mechanical = case.get('mechanics', case)
                joints = mechanical.get('joint_system', {}).get('joints', [])
                structural = mechanical.get('beam_system') is not None or mechanical.get('anchor_system') is not None
                if len(joints) > 1 or structural:
                    targeted.append(dict(path=name, joint_count=len(joints), structural_system=structural))
    files = ['interim/2026-11-16-issue-triage/todo.ja.md', str(FOLDER/'README.ja.md'),
             str(FOLDER/'acceptance.json'), str(FOLDER/'evidence-manifest.json'),
             str(FOLDER/'retained-failures/README.ja.md'), 'src/joints.rs', 'src/stages.rs', 'src/beams.rs']
    sources = {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in files}
    return dict(schema='dry-multi-joint-scope-audit/1', base_commit=PINNED, accepted=None,
                existing_limited_acceptance=acceptance['combination_1'],
                requested_all_combinations_passed=acceptance['requested_all_combinations_passed'],
                archive_sha256=digest, verified_scope_entries=entries, archived_refusals=logs,
                archived_multiple_or_structural_inputs=targeted,
                exact_targeted_archive_input='未確認' if not targeted else 'see list; applicability must be checked',
                verified_source_conditions=[
                    dict(source='src/joints.rs:431-439',
                         condition='any birth_policy AND any crack_tip AND (joints.len()!=1 OR joints[0].water exists)',
                         message='declared-tip JOINT BIRTH is limited to one dry interface',
                         verification='source read only; no new production execution'),
                    dict(source='src/stages.rs:990-998',
                         condition='one JOINT has crack_tip AND birth_policy; beam_system OR anchor_system exists',
                         message='declared-tip JOINT BIRTH with structural members is unverified',
                         verification='source read only; no new production execution')],
                unknown=['old exact dry multiple-JOINT input/log', 'old exact dry beam/plate/anchor composite input/log',
                         'general dry JOINT without crack_tip/BIRTH rejection is NOT established by these guards',
                         'plate-specific combinations not identified from archived refusals'],
                sources_sha256=sources, production_solver_invoked=False,
                A_01='未達', NC_01='未達・対象外', P6='未達・対象外')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument('--archive', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    dump(args.output, audit(args.root, args.archive or args.root/FOLDER/'evidence.zip'))


if __name__ == '__main__':
    main()
