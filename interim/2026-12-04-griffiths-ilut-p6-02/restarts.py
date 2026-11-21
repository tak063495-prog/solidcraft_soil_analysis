"""Real checkpoint files, separate PIDs, full bit comparisons and CLI refusals."""
import argparse
import copy
import json
from pathlib import Path
from study import ROOT, bits, command, fields, read, sha, write


def lines(p):
    with p.open(encoding='utf-8') as f: return [json.loads(x) for x in f]


def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--observer',type=Path,required=True);p.add_argument('--cli',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=a.study.resolve()/'stride5'; inp=source/'input'; obs=a.observer.resolve(); cli=a.cli.resolve()
    prefix=[inp/'mesh.json',inp/'ground.json']; processes={}; traces={};runs={}
    for mode in ['default','ilut']:
        case=inp/('case.json' if mode=='default' else 'case-ilut.json')
        dest=out/mode
        processes[mode]=command([obs,'solve',*prefix,case,dest],out/(mode+'.log'),240)
        if processes[mode]['exit_code']!=0: raise ValueError('trace failed; retained process log')
        traces[mode]=lines(dest/'accepted.jsonl');runs[mode]=read(dest/'run.json')
        assert runs[mode]['complete']
    alignment=[]
    right={(x['stage_index'],x['progress']):x for x in traces['ilut']}
    for x in traces['default']:
        key=(x['stage_index'],x['progress'])
        if key in right:
            alignment.append(dict(stage_index=key[0],progress=key[1],**fields(x['state'],right[key]['state'])))
    restarts=[]
    geo=sum(x['accepted'] for x in runs['ilut']['geostatic']['attempts'])
    checkpoints=[2,geo+2]
    for n in checkpoints:
        cp=out/'ilut'/f'checkpoint-{n}.json';dest=out/f'resume-{n}'
        process=command([obs,'resume',*prefix,inp/'case-ilut.json',dest,cp],out/f'resume-{n}.log',240)
        actual=lines(dest/'accepted.jsonl');expected=[x for x in traces['ilut'] if x['accepted_index']>n]
        restarts.append(dict(checkpoint=str(cp.relative_to(out)),checkpoint_sha256=sha(cp),
            producing_pid=processes['ilut']['pid'],process=process,separate_pid=process['pid']!=processes['ilut']['pid'],
            later_states=len(actual),all_later_states_equal=bits(actual)==bits(expected),
            whole_run_equal=bits(read(dest/'run.json'))==bits(runs['ilut'])))
    # Ordinary CLI issues its own checkpoint; never rebind the diagnostic SHA.
    cp=out/'cli-checkpoint.json'; dest=out/'cli-resumed.json'
    create=command([cli,'checkpoint-stages',*prefix,inp/'case-ilut.json',geo+2,cp],out/'cli-checkpoint.log')
    resume=command([cli,'resume-stages',*prefix,inp/'case-ilut.json',cp,dest],out/'cli-resume.log')
    normal=dict(create=create,resume=resume,separate_pid=create['pid']!=resume['pid'],
                whole_run_equal=bits(read(dest))==bits(read(source/'ilut.json')))
    # Input settings mismatch and invalid values through the ordinary CLI.
    rejections=[]
    def reject(label,case,cp,expected):
        folder=out/label;folder.mkdir();path=folder/'case.json';write(path,case)
        result=folder/'result.json'
        args=[cli,'resume-stages',*prefix,path,cp,result] if cp else [cli,'solve-stages',*prefix,path,result]
        process=command(args,folder/'process.log',60)
        log=(folder/'process.log').read_text()
        row=dict(label=label,process=process,expected=expected,result_exists=result.exists(),
            rejected_as_expected=process['exit_code']==1 and not process['timed_out'] and not result.exists() and expected in log)
        rejections.append(row)
    base=read(inp/'case-ilut.json')
    for label,change,expected in [
        ('drop-negative',{'drop_tolerance':-1.},'ILUT requires'),
        ('drop-too-large',{'drop_tolerance':1.01},'ILUT requires'),
        ('fill-zero',{'fill_limit':0},'ILUT requires'),
        ('fill-too-large',{'fill_limit':4097},'ILUT requires'),
        ('fill-fraction',{'fill_limit':1.5},'invalid type'),
        ('unknown-option',{'ordering':'colamd'},'unknown field')]:
        case=copy.deepcopy(base);case['settings']['linear_solver']['options']['ilut'].update(change)
        reject(label,case,None,expected)
    for label,change in [('mismatch-drop',{'drop_tolerance':1e-3}),('mismatch-fill',{'fill_limit':201}),('mismatch-preconditioner',None)]:
        case=copy.deepcopy(base)
        if change: case['settings']['linear_solver']['options']['ilut'].update(change)
        else: del case['settings']['linear_solver']['options']['ilut']
        reject(label,case,cp,'checkpoint input/solver hash mismatch')
    # Same observer binary, outer checksum recalculated by the Rust canonical hasher.
    saved=read(out/'ilut'/f'checkpoint-{geo+2}.json')
    mutations=[]
    for label in ['checksum','saved-drop','saved-fill','displacement','stress','linear-history']:
        folder=out/('mutate-'+label);folder.mkdir();bad=copy.deepcopy(saved)
        if label in ['saved-drop','saved-fill']:
            option='drop_tolerance' if label=='saved-drop' else 'fill_limit'
            bad['run']['case']['settings']['linear_solver']['options']['ilut'][option]*=2
        elif label=='linear-history':
            bad['run']['geostatic']['attempts'][0]['linear_solves'][0]['iterations']+=1
        elif label=='stress':
            bad['run']['last_committed']['material_states'][0]['effective_stress_compression_positive_kpa'][0]+=1.
        else: bad['run']['last_committed']['displacement_m'][0][0]+=1.
        raw=folder/'raw.json';write(raw,bad);target=raw;rehash=None
        if label!='checksum':
            target=folder/'rehashed.json'
            rehash=command([obs,'rehash',raw,target],folder/'rehash.log')
        dest=folder/'resume'
        process=command([obs,'resume',*prefix,inp/'case-ilut.json',dest,target],folder/'process.log',60)
        expected='checkpoint history/load checksum mismatch' if label=='checksum' else 'checkpoint accepted-path replay mismatch'
        log=(folder/'process.log').read_text()
        mutations.append(dict(label=label,rehash=rehash,process=process,expected=expected,
            rejected_as_expected=process['exit_code']!=0 and not process['timed_out'] and expected in log and not (dest/'run.json').exists()))
    write(out/'restarts.json',dict(processes=processes,restarts=restarts,ordinary_cli=normal,
        accepted_state_counts={k:len(v) for k,v in traces.items()},common_accepted_states=len(alignment),
        default_only=len(traces['default'])-len(alignment),ilut_only=len(traces['ilut'])-len(alignment),
        full_state_field_differences=alignment,input_rejections=rejections,mutations=mutations,
        observer_sha256=sha(obs),cli_sha256=sha(cli),passed=False,p6_complete=False))


if __name__=='__main__': main()
