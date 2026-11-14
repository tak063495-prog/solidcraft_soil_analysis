"""Record explicit refusals; a refusal is never a curved-shell benchmark pass."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import numpy as np
from plate6_offset_study import inputs,write

def probe(output,binary):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    mesh,ground,case,_=inputs();variants=[]
    def add(name,change):
        m,g,c=copy.deepcopy((mesh,ground,case));change(m,g,c);variants.append((name,m,g,c))
    add('nonzero_offset_without_new_formulation',lambda m,g,c:c['beam_system']['plates'][0].update(formulation='mindlin_tri6'))
    add('outside_half_thickness',lambda m,g,c:c['beam_system']['plates'][0].update(attachment_offset_m=.051))
    add('initial_water',lambda m,g,c:g['initialization'].update(water_level_m=.5))
    add('staged_water',lambda m,g,c:c['stages'].append(dict(kind='hydrostatic_water',name='wet',target_level_m=.5)))
    add('joint',lambda m,g,c:c.update(joint_system=dict(schema='soil-fem-joints/1',geometric_tolerance_m=1e-8,maximum_increment_m=.01,tangent_step_m=1e-8,joints=[])))
    add('soil_birth',lambda m,g,c:c['stages'].append(dict(kind='birth',name='new_soil',bodies=['UpperClay'],states={},placement='cad_new_nodes')))
    add('soil_death',lambda m,g,c:c['stages'].append(dict(kind='excavate',name='cut',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.)))
    def plastic_seed(m,g,c):
        for state in g['initialization']['geostatic']['seed']['states'].values():state['equivalent_plastic_strain']=.001
    add('plastic_soil_history',plastic_seed)
    def fold(m,g,c):
        p=c['beam_system']['plates'][1];p['nodes']=[p['nodes'][i] for i in (0,2,1,5,4,3)]
    add('folded_offset_normals',fold)
    for a in range(3):
        def bend(m,g,c,a=a):
            name=c['beam_system']['plates'][0]['nodes'][3]
            node=next(n for n in c['beam_system']['nodes'] if n['id']==name)
            soil=m['node_sets'][node['soil_node_set']][0]
            m['nodes'][soil][a]+=.001;node['position_m']=m['nodes'][soil].copy()
        add('curved_edge_component_'+str(a),bend)
    # Geometric probes on genuine quadratic cylinder/sphere nodal positions.
    # No classical roof/pinch load or mesh-convergence claim is made here.
    rs=np.array([[0.,0.],[.2,0.],[0.,.2],[.1,0.],[.1,.1],[0.,.1]])
    for kind in ('cylinder','sphere'):
        def shell(m,g,c,kind=kind):
            xyz=([[np.sin(r),s,np.cos(r)] for r,s in rs] if kind=='cylinder' else
                 [[np.sin(r)*np.cos(s),np.sin(s),np.cos(r)*np.cos(s)] for r,s in rs])
            cfg=c['beam_system'];cfg['nodes']=[dict(id=str(i),position_m=p,match_soil_node=False) for i,p in enumerate(xyz)]
            cfg['members']=[];cfg['nodal_loads']=[];cfg['rotation_constraints']=[];cfg['translation_constraints']=[]
            cfg['plates']=[dict(id='shell',nodes=list(map(str,range(6))),formulation='mindlin_tri6_offset',attachment_offset_m=0.,section_id='wall_section',initially_active=True)]
            c['stages']=[];c['verification_reject_once_stages']=[]
        add(kind+'_quadratic_geometry',shell)
    records=[]
    for name,m,g,c in variants:
        folder=output/name;folder.mkdir()
        for field,value in [('mesh',m),('ground',g),('stages',c)]:write(folder/(field+'.json'),value)
        p=subprocess.run([str(binary),'solve-stages',*(str(folder/(n+'.json')) for n in ('mesh','ground','stages')),str(folder/'must-not-exist.json')],capture_output=True,text=True)
        assert p.returncode!=0,(name,'unverified scope accepted')
        assert not (folder/'must-not-exist.json').exists()
        records.append(dict(case=name,returncode=p.returncode,stderr=p.stderr,stdout=p.stdout,refused=True))
    result=dict(schema='soil-fem-plate6-curved-offset-rejections/1',records=records,refusal_contract_passed=True,
        curved_tri6_accepted=False,curved_benchmarks_passed=False,
        classical_benchmarks=[dict(name=n,status='not_run_no_validated_curved_kernel',passed=False) for n in ('Scordelis-Lo roof','pinched cylinder','hemispherical shell')])
    write(output/'rejections.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--binary',type=Path,default=Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d')
    args=p.parse_args();print(json.dumps(probe(args.output,args.binary),indent=2))
