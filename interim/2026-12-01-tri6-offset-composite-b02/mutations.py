"""Same producer executable: distinguish checksum, replay and typed refusals."""
import argparse
import copy
from pathlib import Path
from cases import read, write
from study import command, sha


def run(formal, out):
    out.mkdir(parents=True,exist_ok=False)
    binary=next((formal/'solver').iterdir());folder=formal/'beam-positive'
    source=folder/'native/checkpoint-10.json';original=read(source);rows=[]
    base=[folder/(s+'.json') for s in ('mesh','ground','stages')]
    def trial(name,change,rehash=True,expected='checkpoint accepted-path replay mismatch'):
        f=out/name;f.mkdir();bad=copy.deepcopy(original);change(bad)
        write(f/'changed.json',bad)
        if rehash:
            sealing=command([binary,'rehash',f/'changed.json',f/'rehashed.json'],f/'rehash.log',expect=None)
            assert sealing['exit_code']==0
            candidate=f/'rehashed.json'
        else:candidate=f/'changed.json'
        p=command([binary,'resume',*base,f/'resumed',candidate],f/'resume.log',expect=None)
        diagnostic=(f/'resume.log').read_text()
        result=dict(name=name,outer_hash_recomputed=rehash,producer_solver_sha256=sha(binary),checkpoint_sha256=sha(candidate),
                    process=p,diagnostic=diagnostic,expected=expected,run_file_exists=(f/'resumed/run.json').exists(),
                    refused=p['exit_code']!=0 and expected in diagnostic and not (f/'resumed/run.json').exists())
        rows.append(result);write(out/'mutations.json',dict(records=rows,passed=False,b02_complete=False))
    def last(b):return b['run']['last_committed']
    trial('plate-installation-reference',lambda b:last(b)['beams']['plate_reference_global_kinematics'][0].__setitem__(0,.001))
    trial('beam-installation-reference',lambda b:last(b)['beams']['reference_global_kinematics'][0].__setitem__(0,.001))
    trial('soil-displacement',lambda b:last(b)['displacement_m'][8].__setitem__(0,.001))
    trial('structure-rotation',lambda b:last(b)['beams']['node_rotation_rad'][0].__setitem__(1,.001))
    trial('soil-stress',lambda b:last(b)['material_states'][0]['effective_stress_compression_positive_kpa'].__setitem__(0,1.))
    trial('soil-plastic-history',lambda b:last(b)['material_states'][0].__setitem__('equivalent_plastic_strain',.001))
    trial('plate-active',lambda b:last(b)['beams']['active_plates'].__setitem__(0,False))
    trial('soil-active',lambda b:last(b)['active_elements'].__setitem__(0,False))
    trial('applied-force',lambda b:last(b)['applied_nodal_load_kn'].__setitem__(0,.001))
    trial('history-frame-reference',lambda b:b['run']['frames'][1]['state']['beams']['plate_reference_global_kinematics'][0].__setitem__(0,.001))
    trial('known-joint-state',lambda b:last(b).__setitem__('joints',{'points':[]}),expected='checkpoint accepted-path replay mismatch')
    trial('known-anchor-state',lambda b:last(b).__setitem__('anchors',{'anchors':[]}),expected='checkpoint accepted-path replay mismatch')
    trial('known-pore-field',lambda b:last(b).__setitem__('pore_pressure_kpa',[1.]*len(last(b)['material_states'])),expected='checkpoint accepted-path replay mismatch')
    trial('known-water-level',lambda b:last(b).__setitem__('hydrostatic_water_level_m',2.),expected='checkpoint accepted-path replay mismatch')
    trial('known-MATSET-material-ids',lambda b:last(b)['element_material_ids'].__setitem__(0,2),expected='checkpoint accepted-path replay mismatch')
    trial('saved-case-offset',lambda b:b['run']['case']['beam_system']['plates'][0].__setitem__('attachment_offset_m',.04),expected='checkpoint accepted-path replay mismatch')
    trial('saved-case-formulation',lambda b:b['run']['case']['beam_system']['plates'][0].__setitem__('formulation','mindlin_tri6'),expected='checkpoint accepted-path replay mismatch')
    trial('unrehashed-plate-reference',lambda b:last(b)['beams']['plate_reference_global_kinematics'][0].__setitem__(0,.001),False,'checkpoint history/load checksum mismatch')
    # Schema spoofing is deliberately separate from authenticated replay.
    trial('fake-joint-state',lambda b:last(b).__setitem__('joint_offset_state',{}),False,'unknown field')
    trial('fake-water-state',lambda b:last(b).__setitem__('offset_water',{}),False,'unknown field')
    trial('fake-material-state',lambda b:last(b).__setitem__('offset_material','MCC'),False,'unknown field')
    trial('fake-anchor-state',lambda b:last(b).__setitem__('offset_anchor',{}),False,'unknown field')
    result=dict(records=rows,refusal_contract_passed=all(r['refused'] for r in rows),passed=False,b02_complete=False,
                unchanged_checkpoint_sha256=sha(source),same_producer_solver=True)
    write(out/'mutations.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('formal',type=Path);p.add_argument('output',type=Path);a=p.parse_args();run(a.formal.resolve(),a.output.resolve())
