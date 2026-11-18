"""Fresh bounded equilibrium paths; explicit diagnostic entry point only."""
import argparse, copy, hashlib, itertools, json, os
from pathlib import Path
import shutil, subprocess, sys, time
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from a03shear_binding import ROOT, HERE, study, fields, reference, controls, Sink, layout, compare
from a03shear_fixture import fixture, plan, audit, pre_solve_forces

def write(path, value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
def sha(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def command(args, folder, name, minimum=1, reverse=False):
    assert minimum==1 and reverse is False
    start=time.monotonic(); env=dict(os.environ)
    with (folder/(name+'.log')).open('x',encoding='utf-8') as log:
        p=subprocess.Popen(list(map(str,args)),cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        timed=False
        try: code=p.wait(timeout=900)
        except subprocess.TimeoutExpired: timed=True; p.kill(); code=p.wait()
    result=dict(command=list(map(str,args)),pid=p.pid,exit_code=code,timed_out=timed,elapsed_seconds=time.monotonic()-start)
    write(folder/(name+'-process.json'),result)
    if code or timed: raise RuntimeError(('numerically incomplete; retain prefix',result))
    return result

def independent(folder, spec):
    with controls(spec['volume'],False,1):
        model=reference.Model(*[study.read(p) for p in study.paths(folder)])
        write(folder/'independent-layout.json',layout(model))
        with (folder/'independent-fields.jsonl').open('x',encoding='utf-8') as f, (folder/'independent-prepared.jsonl').open('x',encoding='utf-8') as p:
            result=reference.solve(model,spec['n'],Sink(f),Sink(p,work=False))
    for row in result['rows']:
        row.update(newton_iterations=0, direct_linear_solves=1, method='direct linear equilibrium; elastic closed-stick scope')
    result.update(volume_order=spec['volume'],minimum_parts=1,reverse_faces=False)
    write(folder/'independent-result.json',result)

def native(folder, spec, probe, solver):
    n = spec['n']
    proc = command([probe, *study.paths(folder), folder/'native'], folder, 'native', spec['minimum'], spec['reverse'])
    run = study.read(folder/'native/run.json')
    assert run['complete'] and run['provenance']['solver_binary_sha256'] == sha(probe)
    groups = [run['geostatic']['attempts']] + [[a for a in run['attempts'] if a['stage'] == s['name']] for s in run['case']['stages']]
    for attempts in groups:
        assert len(attempts) == n and all(a['accepted'] for a in attempts)
        assert [a['progress_to'] for a in attempts] == [i/n for i in range(1, n+1)]
        assert all(a['relative_force_residual'] <= run['case']['settings']['relative_force_tolerance'] for a in attempts)
        assert all(l['true_relative_residual'] <= 1e-10 for a in attempts for l in a['linear_solves'])
    audits = []
    with controls(5, spec['reverse'], spec['minimum']):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
        # Legacy layout audit writes an independent geometry file before any comparison.
        geometry = study.geometry(model, run, folder)
        (folder/'independent-layout.json').rename(folder/'native-public-layout.json')
        mapping = layout(model)
        write(folder/'native-layout.json', mapping)
        prep = run['transition_audits'][0]
        selected = set(prep['selected_elements'])
        fresh = run['case']['stages'][0]['states']['UpperClay']
        assert all(study.exact(s) == study.exact(fresh) for q, s in zip(run['integration_points'], prep['prepared_state']['material_states']) if q['element_index'] in selected)
        names = ['geostatic'] + [s['name'] for s in run['case']['stages']]
        with (folder/'native-fields.jsonl').open('x', encoding='utf-8') as f:
            sink = Sink(f)
            with (folder/'native/accepted.jsonl').open(encoding='utf-8') as a, (folder/'native/constitutive.jsonl').open(encoding='utf-8') as b:
                count = 0
                for state, diagnostic in itertools.zip_longest(a, b):
                    assert state is not None and diagnostic is not None
                    state, diagnostic = json.loads(state), json.loads(diagnostic)
                    if count == n:
                        sink.placement(prep['prepared_state']['displacement_m'])
                    values, audit = fields.native_fields(model, state, diagnostic, run)
                    assert not np.any(values['material_return']) and not np.any(values['plastic_work'])
                    if count >= 2*n:
                        assert np.all(values['contact_status'] == 1) and not np.any(values['plastic_dissipation'])
                    sink(fields.row(names[count//n], (count % n+1)/n), values)
                    audits.append(dict(row=fields.row(names[count//n], (count % n+1)/n), **audit))
                    count += 1
                assert count == 4*n
        with (folder/'native-prepared.jsonl').open('x', encoding='utf-8') as f:
            sink = Sink(f, work=False)
            for entry in study.read(folder/'native/preparations.json'):
                audit = next(a for a in run['transition_audits'] if a['stage'] == entry['stage'])
                values, check = fields.native_fields(model, audit['prepared_state'], entry['diagnostic'], run, True)
                sink(fields.row(entry['stage'], 0.), values)
    write(folder/'native-audits.json', dict(layout=geometry, rows=audits, fresh_state_exact=True))
    whole = command([probe, *study.paths(folder), folder/'whole', '--whole'], folder, 'whole', spec['minimum'], spec['reverse'])
    assert study.exact(run) == study.exact(study.read(folder/'whole/run.json'))
    standard = None
    if spec['n'] == 8 and spec['minimum'] == 1 and not any(spec[k] for k in ('element', 'face', 'reverse')):
        standard = command([solver, 'solve-stages', *study.paths(folder), folder/'production-run.json'], folder, 'production')
        left, right = copy.deepcopy(run), study.read(folder/'production-run.json')
        for v in (left, right):
            v['provenance'].pop('solver_binary_sha256')
        assert study.exact(left) == study.exact(right)
    return dict(process=proc, whole_process=whole, production_process=standard, states=4*n,
                observed_unobserved_exact_except_elapsed=True, production_default_exact=standard is not None)


def run(out, build):
    out=out.resolve(); out.mkdir(parents=True,exist_ok=False)
    frozen=study.read(HERE/'plan.json'); assert frozen==plan()
    write(out/'plan.json',frozen)
    bindings={}
    for path in sorted(list(HERE.glob('*.py'))+[HERE/'SPEC.ja.md',HERE/'plan.json']):
        bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    write(out/'source-binding.json',bindings)
    shutil.copytree(build/'source',out/'isolated-source')
    shutil.copy2(build/'binding.json',out/'build-binding.json')
    suffix='.exe' if os.name=='nt' else ''
    probe=build/('target/release/a03shear_probe'+suffix); solver=ROOT/('target/release/soil-fem-3d'+suffix)
    write(out/'binaries.json',{k:dict(path=str(p),sha256=sha(p)) for k,p in [('probe',probe),('production',solver)]})
    paths=[]; failures=[]
    # Audit/assemble EVERY fixed input before the first equilibrium solve.
    refusals=[]
    for spec in frozen['specs']:
        folder=out/spec['name']
        try:
            inputs=fixture(spec['n'],spec['height'],spec['stress'])
            study.inputs(folder,inputs); write(folder/'spec.json',spec)
            write(folder/'input-audit.json',audit(*inputs))
            with controls(spec['volume'],False,1):
                write(folder/'pre-solve-seed-forces.json',pre_solve_forces(inputs))
        except (ValueError, AssertionError) as error:
            folder.mkdir(exist_ok=True)
            refusal=dict(path=spec['name'],status='scope_refused_before_solve',error=repr(error),solves=0)
            refusals.append(refusal); write(folder/'scope-refusal.json',refusal)
    write(out/'pre-solve-matrix-audit.json',dict(all_inputs_audited_before_first_solve=True,refusals=refusals))
    for spec in frozen['specs']:
        folder=out/spec['name']
        if any(r['path']==spec['name'] for r in refusals): continue
        for engine in spec['engines']:
            try:
                if engine=='native': value=native(folder,spec,probe,solver)
                else:
                    proc=command([sys.executable,HERE/'a03shear_run.py','--independent',folder],folder,'independent')
                    value=dict(process=proc,states=4*spec['n'])
                paths.append(dict(path=spec['name'],engine=engine,**value))
            except Exception as error:
                failure=dict(path=spec['name'],engine=engine,status='numerically_incomplete',error=repr(error))
                failures.append(failure); write(folder/(engine+'-failure.json'),failure)
        if len(spec['engines'])==2 and not any(f['path']==spec['name'] for f in failures):
            rows=[list(fields.load_fields(folder/(e+'-fields.jsonl'))) for e in spec['engines']]
            write(folder/'implementation-comparison.json',compare(*rows))
        print(spec['name'],'finished',flush=True)
    assert bindings=={p:sha(ROOT/p) for p in bindings}
    write(out/'driver-result.json',dict(passed=False,a03_complete=False,diagnostic_only=True,paths=paths,failures=failures,
        formal_paths=len(paths),formal_states=sum(p['states'] for p in paths),planned_paths=63,scope_refusals=refusals,
        sources_unchanged_after_solves=True,refusals_unchanged_not_resolved=True))
    if failures or refusals: raise RuntimeError('numerically incomplete paths retained')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--build',type=Path);p.add_argument('--independent',type=Path)
    a=p.parse_args()
    if a.independent: independent(a.independent,study.read(a.independent/'spec.json'))
    else: run(a.output,a.build.resolve())
