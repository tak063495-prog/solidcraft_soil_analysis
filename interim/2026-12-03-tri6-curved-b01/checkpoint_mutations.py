"""Authenticate content with the identical producer binary after recomputing outer hash."""
import argparse,copy
from pathlib import Path
from b01_cases import read,write
from study import command,sha


def run(formal,out):
    out.mkdir(parents=True,exist_ok=False);binary=next((formal/'solver').iterdir());rows=[]
    for name in ('attached-n2-s4','roof-n2-s4'):
        folder=formal/name;source=folder/'native/checkpoint-10.json'
        if not source.exists():raise FileNotFoundError(source)
        original=read(source);base=[folder/(k+'.json') for k in ('mesh','ground','stages')]
        def last(x):return x['run']['last_committed']
        def trial(label,fn,rehash=True,expected='checkpoint accepted-path replay mismatch'):
            f=out/name/label;f.mkdir(parents=True);bad=copy.deepcopy(original);fn(bad);write(f/'changed.json',bad)
            if rehash:
                command([binary,'rehash',f/'changed.json',f/'rehashed.json'],f/'rehash.log');candidate=f/'rehashed.json'
            else:candidate=f/'changed.json'
            p=command([binary,'resume',*base,f/'resumed',candidate],f/'resume.log',expect=None)
            log=(f/'resume.log').read_text();row=dict(case=name,name=label,outer_hash_recomputed=rehash,
                solver_sha256=sha(binary),checkpoint_sha256=sha(candidate),process=p,diagnostic=log,expected=expected,
                run_file_exists=(f/'resumed/run.json').exists(),refused=p['exit_code']!=0 and expected in log and not (f/'resumed/run.json').exists())
            rows.append(row);write(out/'mutations.json',dict(passed=False,b01_complete=False,records=rows))
        trial('installation-reference',lambda x:last(x)['beams']['plate_reference_global_kinematics'][0].__setitem__(0,.001))
        trial('rotation',lambda x:last(x)['beams']['node_rotation_rad'][0].__setitem__(1,.001))
        trial('soil-stress',lambda x:last(x)['material_states'][0]['effective_stress_compression_positive_kpa'].__setitem__(0,1.))
        trial('applied-load',lambda x:last(x)['applied_nodal_load_kn'].__setitem__(0,.001))
        trial('history-frame',lambda x:x['run']['frames'][1]['state']['beams']['plate_reference_global_kinematics'][0].__setitem__(0,.001))
        trial('saved-case-e',lambda x:x['run']['case']['beam_system']['plates'][0].__setitem__('attachment_offset_m',.04))
        def normal(x):
            p=x['run']['case']['beam_system']['plates'][0];p['nodes']=[p['nodes'][j] for j in (0,2,1,5,4,3)]
        trial('saved-case-normal',normal)
        def both(x):
            normal(x);p=x['run']['case']['beam_system']['plates'][0];p['attachment_offset_m']=-p.get('attachment_offset_m',0.)
        trial('saved-case-normal-e-pair',both)
        trial('saved-case-curved-midpoint',lambda x:x['run']['case']['beam_system']['nodes'][-1]['position_m'].__setitem__(2,1.234))
        trial('saved-case-outside-half-thickness',lambda x:x['run']['case']['beam_system']['plates'][0].__setitem__('attachment_offset_m',.501))
        trial('checksum-not-updated',lambda x:last(x)['beams']['node_rotation_rad'][0].__setitem__(1,.001),False,'checkpoint history/load checksum mismatch')
        trial('fake-per-node-normal',lambda x:last(x)['beams'].__setitem__('node_normals',[[0,0,1]]),False,'unknown field')
    result=dict(passed=False,b01_complete=False,records=rows,refusal_contract_passed=all(r['refused'] for r in rows),same_producer_solver=True)
    write(out/'mutations.json',result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
