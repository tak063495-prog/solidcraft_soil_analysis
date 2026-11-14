"""Prepare an explicitly identified user section; never invent geometry or run SRM.

Uses the unchanged durable production job and the read-only Rust Fs verifier.
Only an actual completed, independently audited search may pass verification.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_new(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def prepare(solver, mesh, ground, case, spec, source, section_id, output):
    if not section_id.strip():
        raise ValueError('実断面を特定するsection-idが必要です')
    inputs = [Path(p).resolve(strict=True) for p in (solver, mesh, ground, case, spec, source)]
    if not all(p.is_file() for p in inputs):
        raise ValueError('全入力は実在するファイルである必要があります')
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    # The source document records the user's provenance assertion; it is not a
    # certification that the geometry matches a surveyed section.
    shutil.copyfile(inputs[5], output / 'section-source')
    subprocess.run([str(inputs[0]), 'prepare-srm-search-job', *map(str, inputs[1:5]),
                    '1', str(output / 'job')], check=True, capture_output=True, text=True)
    manifest = {
        'schema': 'soil-fem-real-section-srm-preparation/1',
        'section_id': section_id, 'geometry_provenance': 'user_supplied_document_not_independently_surveyed',
        'source_document_sha256': digest(output / 'section-source'),
        'job_manifest_sha256': digest(output / 'job' / 'job.json'),
        'search_executed_by_preparer': False, 'exact_fs_for_this_case': False,
        'real_section_srm_complete': False,
    }
    write_new(output / 'preparation.json', manifest)
    return manifest


def verify(bundle, output):
    bundle = Path(bundle).resolve(strict=True)
    m = json.loads((bundle / 'preparation.json').read_text(encoding='utf-8'))
    if (m['schema'] != 'soil-fem-real-section-srm-preparation/1'
            or not m['section_id'].strip()
            or m['source_document_sha256'] != digest(bundle / 'section-source')
            or m['job_manifest_sha256'] != digest(bundle / 'job' / 'job.json')
            or m['search_executed_by_preparer'] is not False
            or m['exact_fs_for_this_case'] is not False
            or m['real_section_srm_complete'] is not False):
        raise ValueError('準備記録・出典・ジョブ結合が変更されています')
    job = bundle / 'job'
    solver = job / 'solver' / 'soil-fem-3d.exe'
    # Even on Unix the durable job's producer binary retains this filename.
    subprocess.run([str(solver), 'verify-srm-search-job', str(job)], check=True,
                   capture_output=True, text=True)
    status = json.loads((job / 'status.json').read_text(encoding='utf-8'))
    if (status['state'] != 'completed' or status['report'] != 'report.json'
            or status['report_sha256'] != digest(job / 'report.json')):
        raise ValueError('探索未完了、または保存時の探索結果SHA-256と不一致です')
    # Missing report, incomplete jobs and altered inputs cannot pass. Exit 2
    # from this verifier denotes a valid but unresolved result, never success.
    subprocess.run([str(solver), 'verify-exact-fs', str(job / 'inputs' / 'mesh.json'),
                    str(job / 'inputs' / 'ground.json'), str(job / 'inputs' / 'case.json'),
                    str(job / 'inputs' / 'specification.json'), str(job / 'report.json'), status['report_sha256'],
                    str(Path(output).resolve())], check=True, capture_output=True, text=True)
    return json.loads(Path(output).read_text(encoding='utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    for arg in ('solver', 'mesh', 'ground', 'case', 'spec', 'source', 'section-id', 'output'):
        p.add_argument('--' + arg, required=True)
    p = sub.add_parser('verify')
    p.add_argument('bundle')
    p.add_argument('output')
    a = parser.parse_args()
    try:
        if a.command == 'prepare':
            result = prepare(a.solver, a.mesh, a.ground, a.case, a.spec, a.source, a.section_id, a.output)
        else:
            result = verify(a.bundle, a.output)
        print(json.dumps(result, ensure_ascii=False, indent=2))
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'未完了または入力・証拠不正: {exc}\n')


if __name__ == '__main__':
    main()
