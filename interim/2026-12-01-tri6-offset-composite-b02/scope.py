"""Actual CLI scope/type/geometry probes, including e=0 guard disguises."""
import argparse
import copy
from pathlib import Path
from cases import ROOT, inputs, read, write
from study import command, sha
import fem_reference as ref
import numpy as np


def run(out, cli):
    out.mkdir(parents=True, exist_ok=False); rows=[]
    m,g,c,_=inputs()
    xyz=np.asarray(m['nodes']);topology=ref.topology(m)
    m['face_sets']['bottom']=[list(face) for face,owners in topology.items() if len(owners)==1 and np.all(xyz[list(face),2]==0.)]
    def check(name, mm, gg, cc, expected, entry='solve-stages', category='scope'):
        f=out/name;f.mkdir()
        for key,v in [('mesh',mm),('ground',gg),('case',cc)]:write(f/(key+'.json'),v)
        target=f/'result.json'
        p=command([cli,entry,f/'mesh.json',f/'ground.json',f/'case.json',target],f/'cli.log',expect=None)
        diagnostic=(f/'cli.log').read_text()
        record=dict(name=name,entry=entry,category=category,expected=expected,diagnostic=diagnostic,process=p,
                    result_file_exists=target.exists(),matched_expected=expected in diagnostic,
                    input_sha256={k:sha(f/(k+'.json')) for k in ('mesh','ground','case')})
        record['refused']=p['exit_code']==1 and expected in diagnostic and not target.exists()
        if target.exists():
            r=read(target);record['complete']=r.get('complete');record['failure_reason']=r.get('failure_reason')
            record['incomplete_scope_response']=p['exit_code']==1 and expected in diagnostic and not r.get('complete',True)
        rows.append(record);write(out/'refusals.json',dict(passed=False,b02_complete=False,records=rows))
        return record
    # Real nonempty JOINT configuration on its own split mesh. Add an actual
    # offset plate independently, leaving the copied JOINT declaration intact.
    jm=read(ROOT/'examples/joint-column-tet10-mesh.json')
    jc=read(ROOT/'examples/stages-joint-verification.json')
    jg=read(ROOT/'examples/ground-joint-mc.json')
    # JOINT uses its own conforming split Tet10 mesh and original nonempty
    # declaration; the offset plate is a real six-node exterior trace. The
    # early scope refusal still does not validate hypothetical future physics.
    for e,label in [(.05,'nonzero'),(0.,'zero')]:
        for what in ('joint','anchor','hydrostatic-initial','hydrostatic-stage','steady-darcy','soil-death','soil-birth','soil-matset','MC-plastic-history','MCC','IDH'):
            mm,gg,cc=copy.deepcopy((m,g,c))
            for p in cc['beam_system']['plates']:p['attachment_offset_m']=e
            expected='offset Tri6 scope excludes'
            if what=='joint':
                mm=copy.deepcopy(jm);gg=copy.deepcopy(jg);cc=copy.deepcopy(jc)
                cc['stages']=[];cc['verification_reject_once_stages']=[]
                topology=ref.topology(mm);face=mm['face_sets']['exterior'][0]
                ids=ref.face_nodes(mm,topology[tuple(sorted(face))][0],face).tolist()
                cfg=copy.deepcopy(c['beam_system']);cfg['members']=[];cfg['nodal_loads']=[]
                cfg['nodes']=[dict(id=f'joint-offset-{n}',position_m=mm['nodes'][n],soil_node_set=f'joint-offset:{n}') for n in ids]
                for n in ids:mm['node_sets'][f'joint-offset:{n}']=[n]
                plate=copy.deepcopy(cfg['plates'][0]);plate.update(id='joint-offset-face',nodes=[n['id'] for n in cfg['nodes']],initially_active=True,attachment_offset_m=e)
                cfg['plates']=[plate];cc['beam_system']=cfg
            elif what=='anchor':
                cfg=read(ROOT/'examples/stages-anchor-slope.json')['anchor_system'];a=cfg['anchors'][0]
                a['head']['position_m']=[.1,.1,.9]
                a['bonded_path']=[{'position_m':[.3,.3,.7]},{'position_m':[.5,.5,.5]}]
                cc['anchor_system']=cfg
            elif what=='hydrostatic-initial':gg['initialization']['water_level_m']=2.;expected='offset Tri6 with water is unverified'
            elif what=='hydrostatic-stage':cc['stages'].append(dict(kind='hydrostatic_water',name='wet',target_level_m=2.))
            elif what=='steady-darcy':
                cc['stages'].append(dict(kind='steady_seepage',name='Darcy',flow=read(ROOT/'examples/flow-seepage-column.json'),exterior_water_level_m=2.))
            elif what=='soil-death':cc['stages'].append(dict(kind='excavate',name='soil-death',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.))
            elif what=='soil-birth':
                cc['initially_inactive_bodies']=['UpperClay']
                cc['stages'].insert(0,dict(kind='birth',name='soil-birth',bodies=['UpperClay'],states={'UpperClay':copy.deepcopy(gg['initialization']['geostatic']['seed']['states']['UpperClay'])},placement='cad_new_nodes'))
            elif what=='soil-matset':cc['stages'].append(dict(kind='material_set',name='soil-matset',bodies=['UpperClay'],material_id=1,history={'mode':'preserve'}))
            elif what=='MC-plastic-history':
                for seed in gg['initialization']['geostatic']['seed']['states'].values():seed['equivalent_plastic_strain']=.001
                expected='offset Tri6 scope requires elastic MC soil'
            elif what=='MCC':
                gg=read(ROOT/'examples/ground-water-column-mcc.json');gg['initialization']['water_level_m']=None
                expected='offset Tri6 soil scope is dry elastic MC only'
            elif what=='IDH':gg=read(ROOT/'examples/ground-hardening-column.json');expected='offset Tri6 soil scope is dry elastic MC only'
            check(label+'-'+what,mm,gg,cc,expected)
        # Existing transient ONE_WAY controller, no invented adapter.
        ow=read(ROOT/'examples/one-way-mcc.json');ow['mechanics']=copy.deepcopy(c)
        for p in ow['mechanics']['beam_system']['plates']:p['attachment_offset_m']=e
        # Existing validation order rejects structures before any transfer.
        mm=copy.deepcopy(m)
        topology=ref.topology(mm);xyz=np.asarray(mm['nodes'])
        mm['face_sets']['bottom']=[list(face) for face,owners in topology.items() if len(owners)==1 and np.all(xyz[list(face),2]==0.)]
        check(label+'-transient-Richards',mm,g,ow,'ONE_WAY requires dry-reference weights',entry='solve-one-way')
    variants=[]
    def add(name,change,expected,category='geometry'):
        mm,gg,cc=copy.deepcopy((m,g,c));change(mm,gg,cc);variants.append((name,mm,gg,cc,expected,category))
    add('outside-positive',lambda m,g,c:c['beam_system']['plates'][0].update(attachment_offset_m=.051),'offset')
    add('outside-negative',lambda m,g,c:c['beam_system']['plates'][0].update(attachment_offset_m=-.051),'offset')
    add('old-formulation-disguise',lambda m,g,c:c['beam_system']['plates'][0].update(formulation='mindlin_tri6'),'offset')
    add('old-alias-disguise',lambda m,g,c:c['beam_system']['plates'][0].update(thickness_offset_m=.05),'unknown field','typed')
    def normal(m,g,c):
        p=c['beam_system']['plates'][1];p['nodes']=[p['nodes'][j] for j in (0,2,1,5,4,3)]
    add('inconsistent-normal',normal,'offset')
    add('inconsistent-sign',lambda m,g,c:c['beam_system']['plates'][1].update(attachment_offset_m=-.05),'offset')
    def curve(m,g,c):
        n=c['beam_system']['nodes'][3];i=m['node_sets'][n['soil_node_set']][0]
        m['nodes'][i][2]+=.001;n['position_m']=m['nodes'][i].copy()
    add('declared-flat-actual-curved',curve,'Mindlin Tri6')
    def mismatch(m,g,c):c['beam_system']['nodes'][0]['position_m'][2]+=.001
    add('declared-attachment-actual-mismatch',mismatch,'Mindlin Tri6')
    add('fake-joint-field',lambda m,g,c:c.update(joints=[]),'unknown field','typed')
    add('fake-water-field',lambda m,g,c:c.update(water_level_m=2.),'unknown field','typed')
    add('fake-anchor-field',lambda m,g,c:c.update(anchors=[]),'unknown field','typed')
    add('fake-material-field',lambda m,g,c:c.update(soil_model='MCC'),'unknown field','typed')
    add('fake-birth-field',lambda m,g,c:c.update(birth_bodies=['UpperClay']),'unknown field','typed')
    add('fake-death-field',lambda m,g,c:c.update(death_bodies=['UpperClay']),'unknown field','typed')
    add('fake-matset-field',lambda m,g,c:c.update(material_set={'UpperClay':1}),'unknown field','typed')
    for name,mm,gg,cc,expected,category in variants:check(name,mm,gg,cc,expected,category=category)
    # A newly entered plastic branch must fail even with no initial plastic seed.
    mm,gg,cc=copy.deepcopy((m,g,c))
    for material in gg['materials']:material['mohr_coulomb']['cohesion_kpa']=.001
    response=check('MC-new-plastic-response',mm,gg,cc,'geostatic cutback exhausted',category='numerical_incompletion_not_scope_refusal')
    response['refused']=False
    response['diagnostic_scope_refusal_confirmed']=False
    response['incomplete_scope_response']=False
    response['judgment']='未達。数値未完了を専用の塑性応答拒否成功へ数えない'
    result=dict(passed=False,b02_complete=False,records=rows,
                refusal_contract_passed=all(r['refused'] for r in rows if r['category']!='numerical_incompletion_not_scope_refusal'),
                numerical_incomplete_cases=[r['name'] for r in rows if r['category']=='numerical_incompletion_not_scope_refusal'],
                note='Early scope rejection does not validate hypothetical combined ownership, initialisation or future physics. Typed fake fields are separate. A written incomplete Run is not a success.')
    write(out/'refusals.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--cli',type=Path,required=True);a=p.parse_args();run(a.output.resolve(),a.cli.resolve())
