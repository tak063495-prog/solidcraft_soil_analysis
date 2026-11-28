"""Explicit opt-in small A-03 diagnostics. Default solvers never import this."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):os.environ[key]='1'
from a03cause_binding import ROOT, OLD, study, fields
from a03cause_reference import independent


def write(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False);f.write('\n')


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f, 'sha256').hexdigest()


def reverse_keys(value):
    if isinstance(value, dict):return {k:reverse_keys(v) for k,v in reversed(list(value.items()))}
    if isinstance(value, (list, tuple)):return [reverse_keys(v) for v in value]
    return value


def source_binding(out):
    tracked=subprocess.check_output(['git','ls-files','src','tools','tests','examples','Cargo.toml','Cargo.lock','rust-toolchain.toml','.github'],cwd=ROOT,text=True).splitlines()
    prior=[p.relative_to(ROOT).as_posix() for p in OLD.glob('*.py')]
    own=[p.relative_to(ROOT).as_posix() for p in Path(__file__).parent.iterdir() if p.suffix=='.py' or p.name=='SPEC.ja.md']
    paths=sorted(set(tracked+prior+own+['tests/test_a03_wet_soil_birth_cause_diagnosis.py']))
    binding={p:sha(ROOT/p) for p in paths};write(out/'source-binding.json',binding)
    archive=out/'source-snapshot';archive.mkdir()
    for p in paths:
        target=archive/p;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,target)
    return binding


def run(out, counts):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    binding=source_binding(out);binpath=out/'bin';binpath.mkdir();binaries={};suffix='.exe' if os.name=='nt' else ''
    for name,source in [('solver',ROOT/('target/release/soil-fem-3d'+suffix)),('probe',ROOT/('target/release/examples/wet_birth_a03_probe'+suffix)),('rehash',ROOT/('target/release/examples/joint_ext_checkpoint_rehash'+suffix))]:
        dest=binpath/(name+suffix);shutil.copy2(source,dest);binaries[name]=dict(path=str(dest),sha256=sha(dest))
    write(out/'binaries.json',binaries)
    probe,solver,rehash=(binaries[n]['path'] for n in ('probe','solver','rehash'))
    report=dict(passed=False,a03_complete=False,diagnostic_only=True,cases={},refusals=[],fresh_paths=[])
    def independent_process(folder,n,order=5,reverse=False):
        command=[sys.executable,Path(__file__),'--independent',folder,'--steps',str(n),'--volume-order',str(order)]
        if reverse:command.append('--reverse-volume')
        proc=study.invoke(command,folder/'independent.log',timeout=600);write(folder/'independent-process.json',proc)
        if proc['exit_code']!=0 or proc['timed_out']:raise RuntimeError(('independent incomplete',folder,proc))
        report['fresh_paths'].append(dict(kind='input-only-independent',folder=folder.name,states=4*n,volume_order=order,reversed_volume_accumulation=reverse,process=proc))
    for n in counts:
        for label,level in [('h1_5',1.5),('h2_0',2.)]:
            name=f'{label}-n{n}';folder=out/name;study.inputs(folder,study.fixture(n,level))
            observed=study.native(folder,Path(probe),Path(solver),n)
            for kind,filename in [('observed-native','native-process.json'),('unobserved-native','whole-process.json'),('production-cli','production-process.json')]:
                report['fresh_paths'].append(dict(kind=kind,folder=name,states=4*n,process=study.read(folder/filename)))
            independent_process(folder,n)
            comparison=fields.compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl')
            prepared=fields.compare(folder/'native-prepared.jsonl',folder/'independent-prepared.jsonl')
            write(folder/'legacy-comparison.json',comparison);write(folder/'legacy-prepared-comparison.json',prepared)
            for key in ('u','stress','traction'):
                assert comparison['all_fields'][key]['relative_l2']<1e-7,(name,key,comparison['all_fields'][key])
            for key in ('internal','applied_load','effective_internal','pore_internal','contact_internal','gravity_load','surcharge_load','exterior_water_load','joint_water_load','holding_load'):
                assert comparison['all_fields'][key]['maximum_absolute']<1e-10,(name,key)
            assert all(v['mismatches']==0 for v in comparison['all_fields'].values() if 'mismatches' in v)
            report['cases'][name]=dict(native=observed,comparison=comparison,preparation=prepared,independent=study.read(folder/'independent-result.json'))
            if n==4:
                report['cases'][name]['restart']=study.restart(folder,Path(probe),n)
                report['cases'][name]['tampering']=study.tamper(folder,Path(probe),Path(rehash),n)
                quad=out/(name+'-gl7');study.inputs(quad,study.fixture(n,level));independent_process(quad,n,7)
                reverse=out/(name+'-reverse-volume');study.inputs(reverse,study.fixture(n,level));independent_process(reverse,n,5,True)
                ordered=out/(name+'-key-order');study.inputs(ordered,reverse_keys(study.fixture(n,level)))
                independent_process(ordered,n)
                proc=study.invoke([probe,*study.paths(ordered),ordered/'whole','--whole'],ordered/'whole.log',timeout=600);write(ordered/'whole-process.json',proc)
                assert proc['exit_code']==0 and not proc['timed_out']
                report['fresh_paths'].append(dict(kind='key-order-native',folder=ordered.name,states=4*n,process=proc))
                original=study.read(folder/'native/run.json');reordered=study.read(ordered/'whole/run.json')
                report['cases'][name]['key_order_full_run_exact_except_elapsed']=study.exact(original)==study.exact(reordered)
                assert study.exact(original['transition_audits'])==study.exact(reordered['transition_audits'])
                assert study.exact(original['last_committed'])==study.exact(reordered['last_committed'])
                assert list(fields.load_fields(folder/'independent-fields.jsonl'))[0][0]==list(fields.load_fields(ordered/'independent-fields.jsonl'))[0][0]
    if 4 in counts:
        report['refusals']=study.refusals(out,Path(solver),4)
        assert len(report['refusals'])==24
        # A surface-crossing element can have every interior GP below water.
        # Do not replace the all-node scope guard by a GP-only saturation test.
        m,g,c=study.fixture(4)
        level=math.nextafter(max(p[2] for p in m['nodes']),-math.inf)
        g['initialization']['water_level_m']=level;c['stages'][0]['hydraulic_initialization']['level_m']=level
        edge=out/'partial-node-above-all-GP-below';study.inputs(edge,(m,g,c))
        proc=study.invoke([solver,'solve-stages',*study.paths(edge),edge/'forbidden.json'],edge/'cli.log',timeout=600)
        write(edge/'process.json',proc)
        assert proc['exit_code']==1 and not proc['timed_out'] and not (edge/'forbidden.json').exists()
        assert 'wet BIRTH requires one fully submerged fresh soil body' in (edge/'cli.log').read_text()
        positions=[p[2] for e in range(len(m['tetrahedra'])) for p in study.fem.volume(m,e)['position']]
        assert max(positions)<level
        report['additional_refusals']=[dict(name=edge.name,explicit_matching_level_m=level,process=proc,result_file_exists=False,
            volume_GP_count=len(positions),maximum_volume_GP_z_m=max(positions),mesh_nodes_above_level=sum(p[2]>level for p in m['nodes']),
            reason='all interior GP below level does not imply all selected nodes submerged; existing scope guard retained')]
    report['fixed_sources_unchanged_after_solves']=binding=={p:sha(ROOT/p) for p in binding}
    assert report['fixed_sources_unchanged_after_solves']
    write(out/'driver-result.json',report)
    print(json.dumps(dict(output=str(out),fresh_paths=len(report['fresh_paths']),refusals=len(report['refusals']),a03_complete=False)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path);p.add_argument('--counts',type=int,nargs='+',default=[2,4,8])
    p.add_argument('--independent',type=Path);p.add_argument('--steps',type=int,default=4);p.add_argument('--volume-order',type=int,choices=(5,7),default=5)
    p.add_argument('--reverse-volume',action='store_true')
    a=p.parse_args()
    if a.independent:independent(a.independent,a.steps,a.volume_order,a.reverse_volume);return
    if a.output is None or any(n<2 or n%2 for n in a.counts) or len(set(a.counts))!=len(a.counts):p.error('unique positive even counts and a new output are required')
    run(a.output,a.counts)


if __name__=='__main__':main()
