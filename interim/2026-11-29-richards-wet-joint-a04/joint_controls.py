"""Fresh existing dry/hydrostatic controls. Each worker has isolated imports."""
import argparse
import copy
import importlib.util
from pathlib import Path
import shutil
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from study import read,write,sha,invoke,exact


def control(out,kind,steps):
    old=ROOT/'interim'/('2026-11-25-joint-birth-projection-a02' if kind=='dry' else '2026-11-26-wet-soil-birth-darcy-a03')
    sys.path.insert(0,str(old));spec=importlib.util.spec_from_file_location('a04_control_'+kind,old/'study.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
    out.mkdir(parents=True,exist_ok=False);binpath=out/'bin';binpath.mkdir();suffix='.exe' if sys.platform=='win32' else ''
    original=ROOT/'target/release/examples'/('joint_a02_history_probe'+suffix if kind=='dry' else 'wet_birth_a03_probe'+suffix)
    probe=binpath/('probe'+suffix);solver=binpath/('solver'+suffix);rehash=binpath/('rehash'+suffix)
    for a,b in ((original,probe),(ROOT/('target/release/soil-fem-3d'+suffix),solver),(ROOT/('target/release/examples/joint_ext_checkpoint_rehash'+suffix),rehash)):shutil.copy2(a,b)
    write(out/'binaries.json',{p.name:sha(p) for p in (probe,solver,rehash)})
    folder=out/'case'
    if kind=='dry':
        s.make_inputs(folder,s.fixture('stress_free',steps,True));s.native(folder,probe,steps)
        command=[sys.executable,old/'study.py','--output',out,'--steps',steps,'--independent',folder]
    else:
        s.inputs(folder,s.fixture(steps,2.));s.native(folder,probe,solver,steps)
        command=[sys.executable,old/'study.py','--output',out,'--steps',steps,'--independent',folder]
    proc=invoke(command,folder/'independent-worker.log');assert proc['exit_code']==0 and not proc['timed_out']
    comparison=s.compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl');write(folder/'comparison.json',comparison)
    for key in ('u','stress','traction'):assert comparison['all_fields'][key]['relative_l2']<1e-7
    for key in ('internal','applied_load'):assert comparison['all_fields'][key]['maximum_absolute']<1e-10
    assert all(v.get('mismatches',0)==0 for v in comparison['all_fields'].values())
    restarts=s.restart(folder,probe,steps);tampering=s.tamper(folder,probe,rehash,steps)
    # Add explicit open/closed and generation spoofing at an accepted checkpoint.
    cps=sorted((folder/'native').glob('checkpoint-*.json'));cp=cps[len(cps)//2];base=read(cp);extra=out/'state-spoof';extra.mkdir();spoof=[]
    for mode in ('open-status','generation','physical-time'):
        bad=copy.deepcopy(base);state=bad['run']['last_committed']
        if mode=='open-status':state['joints']['points'][0]['status']='open' if state['joints']['points'][0]['status']!='open' else 'stick'
        elif mode=='generation':state['joints']['lifecycle'][0]['generation']+=1
        else:state['time_s']=1.25
        raw=extra/(mode+'-raw.json');sealed=extra/(mode+'-sealed.json');write(raw,bad)
        seal=invoke([rehash,raw,sealed],extra/(mode+'-rehash.log'))
        if mode=='physical-time':
            assert seal['exit_code']==1 and not seal['timed_out'] and not sealed.exists();assert 'unknown field' in (extra/(mode+'-rehash.log')).read_text()
            spoof.append(dict(mode=mode,reason='quasistatic state has no physical time; typed input rejects unknown field',process=seal,outer_hash_recomputed=False));continue
        assert seal['exit_code']==0
        result=invoke([probe,folder/'mesh.json',folder/'ground.json',folder/'case.json',extra/(mode+'-resume'),'--resume',sealed],extra/(mode+'.log'))
        text=(extra/(mode+'.log')).read_text();assert result['exit_code']==1 and not result['timed_out'] and 'checkpoint accepted-path replay mismatch' in text,text
        spoof.append(dict(mode=mode,reason=text,process=result,outer_hash_recomputed=True,producer_binary_sha256=sha(probe)))
    prepared=None
    if kind=='hydrostatic':
        prepared=s.compare(folder/'native-prepared.jsonl',folder/'independent-prepared.jsonl');write(folder/'prepared-comparison.json',prepared)
    result=dict(passed=False,a04_complete=False,kind=kind,judgment='既存の限定受入の維持・再確認',
        scope='A-02 dry single tip B/death/rebirth' if kind=='dry' else 'A-03 fully submerged fresh MC/CAD matching hydrostatic head/CONTINUOUS B',
        fresh_native_paths=1,fresh_input_only_independent_paths=1,comparison=comparison,prepared=prepared,restarts=restarts,tampering=tampering,state_spoofing=spoof,
        hydraulic_storage_m3=None,hydraulic_flow_m3_s=None,mass_balance=None,physical_time_s=None,
        not_a_Richards_joint_solve=True,independent_process=proc)
    write(out/'measurements.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--kind',choices=('dry','hydrostatic'),required=True);p.add_argument('--steps',type=int,default=4);a=p.parse_args();control(a.output.resolve(),a.kind,a.steps)
