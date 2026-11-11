"""Audit every saved material state and native GP enum; report switching brackets."""
import argparse
from pathlib import Path
import json
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import mapping as m
c=m.c

def run(out,target):
    materials={};switches={}
    for folder in sorted(out.glob('h*-g*')):
        raw=folder/'native-raw/accepted.jsonl'
        if not (folder/'native/manifest.json').exists():continue
        h=int(folder.name.split('-')[0][1:]);g=int(folder.name.split('-')[1][1:]);mesh=c.read(out/f'mesh-h{h}.json')
        nmat=27*len(mesh['tetrahedra']);ngp=32*4**h;count=0;prev=None;events=[]
        with raw.open() as stream:
            for i,line in enumerate(stream):
                state=json.loads(line);ms=state['material_states'];assert len(ms)==nmat
                for x in ms:
                    assert x['equivalent_plastic_strain']==x['plastic_work_kpa']==0
                    assert not any(x['plastic_strain_compression_positive']) and x['last_return']=='elastic'
                    assert set(x)=={'effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','equivalent_plastic_strain','plastic_work_kpa','last_return'}
                js=state['joints']['points'];assert len(js)==ngp
                status=[x['status'] for x in js];gap=np.array([np.dot(x['jump_m'],x['normal']) for x in js])
                assert all(s in ('open','stick','slip') for s in status)
                if prev is not None:
                    oldgap,oldstatus,oldprogress,oldstage=prev
                    phase,step=divmod(i,g);progress=(step+1)/g;stage=c.STAGES[phase]
                    for q in np.flatnonzero((gap>0)!=(oldgap>0)):
                        events.append(dict(point=int(q),stage=stage,progress_from=oldprogress if oldstage==stage else 0.,progress_to=progress,old_gap_m=float(oldgap[q]),new_gap_m=float(gap[q]),old_status=oldstatus[q],new_status=status[q],kind='opening' if gap[q]>0 else 'closing'))
                phase,step=divmod(i,g);prev=(gap,status,(step+1)/g,c.STAGES[phase]);count+=1
        assert count==7*g
        materials[folder.name]=dict(material_states_checked=count*nmat,all_elastic_return=True,equivalent_plastic_strain_zero=True,plastic_work_zero=True,plastic_strain_zero=True,native_contact_gp_states_checked=count*ngp)
        s=c.old.Water(mesh,c.read(out/'ground.json'),c.read(out/'stages.json')).contact
        switches[folder.name]=dict(points=s.points,events=events,strict_gap_sign_definition=True,switch_times_are_brackets_not_exact=True)
    maps=m.Maps([c.read(out/f'mesh-h{h}.json') for h in range(3)],c.read(out/'stages.json'))
    common={}
    for folder in sorted(out.glob('h*-g*')):
        if not (folder/'native/manifest.json').exists():continue
        h=int(folder.name.split('-')[0][1:]);g=int(folder.name.split('-')[1][1:]);brackets=[None]*len(maps.contact_position)
        prev=None;prevrow=None
        for row,f in c.old.states(folder/'native'):
            gap=maps.operators[h]['contact']@f['gap']
            if row['stage']=='open':
                for q in np.flatnonzero(gap>0):
                    if brackets[q] is None:
                        brackets[q]=dict(progress_from=prevrow['progress'] if prevrow['stage']=='open' else 0.,progress_to=row['progress'],gap_from_m=float(prev[q]),gap_to_m=float(gap[q]))
            prev=gap;prevrow=row
        common[folder.name]=brackets
    comparisons={}
    for g in (32,64,256):
        for h in (0,1):
            names=[f'h{h}-g{g}',f'h{h+1}-g{g}']
            if not all(n in common for n in names):continue
            pairs=[(a,b) for a,b in zip(*(common[n] for n in names)) if a is not None and b is not None]
            delta=[.5*(a['progress_from']+a['progress_to']-b['progress_from']-b['progress_to']) for a,b in pairs]
            comparisons['__'.join(names)]=dict(points_with_both_openings=len(pairs),maximum_midpoint_progress_difference=float(np.max(abs(np.array(delta)),initial=0)),nonoverlapping_brackets=sum(max(a['progress_from'],b['progress_from'])>min(a['progress_to'],b['progress_to']) for a,b in pairs),exact_switch_time=False)
    c.write(target/'common-opening-brackets.json',dict(position_m=maps.contact_position,root_face=maps.contact_root,first_opening_in_open_stage=common,comparisons=comparisons,strict_gap_sign=True,time_alignment_used=False))
    c.write(target/'raw-material-audit.json',materials);c.write(target/'switch-brackets.json',switches)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.input,a.output)
