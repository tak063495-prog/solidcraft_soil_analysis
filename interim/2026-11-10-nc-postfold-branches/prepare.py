"""Verify prior inputs and build the existing observer from this exact source tree."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def run(out):
    out.mkdir(parents=True,exist_ok=False)
    start=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    if subprocess.run(['git','merge-base','--is-ancestor','b235824ea55a0e70b8788e586ee8562549f8baf5','HEAD'],cwd=ROOT).returncode:
        raise ValueError('specified PR tip is not an ancestor')
    previous=ROOT/'interim/2026-11-10-nc-tangent-branch-diagnosis'
    extraction=module('branch_previous_extract',previous/'extract.py').run(out/'prior')
    module('branch_previous_observer',previous/'prepare_observer.py').run(out/'observer')
    original=json.loads((out/'observer/preparation.json').read_text(encoding='utf-8'))
    for name,sha in original['source_original_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=sha:
            raise ValueError('current source differs from observer input')
    provenance=dict(invocation_commit=start,required_ancestor='b235824ea55a0e70b8788e586ee8562549f8baf5',
        prior_archive_sha256=extraction['archive_sha256'],observer_binary_sha256=original['binary_sha256'],
        observer_source_matches_current=True,historical_helper_commit_label_is_not_invocation_commit=original['starting_commit'],
        production_modified=False,numerical_acceptance=False)
    (out/'preparation.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
