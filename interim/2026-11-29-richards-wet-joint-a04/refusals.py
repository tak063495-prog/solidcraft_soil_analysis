"""Actual CLI fail-closed scope probes, never a new hydraulic/JOINT adapter."""
import argparse
import copy
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from study import ROOT,fixture,read,write,invoke


def run(out,solver):
    out.mkdir(parents=True,exist_ok=False);rows=[]
    def check(name,command,mesh,ground,case,expected):
        f=out/name;f.mkdir()
        for k,v in zip(('mesh','ground','case'),(mesh,ground,case)):write(f/(k+'.json'),v)
        target=f/'result.json';cmd=[solver,command,f/'mesh.json']+([f/'ground.json'] if command!='solve-richards' else [])+[f/'case.json',target]
        p=invoke(cmd,f/'cli.log');text=(f/'cli.log').read_text()
        assert p['exit_code']==1 and not p['timed_out'] and not target.exists(),(name,p,text)
        assert expected in text,(name,expected,text)
        rows.append(dict(name=name,expected=expected,diagnostic=text,process=p,result_file_exists=False,judgment='未達／未対応（現行拒否を保持）'))
    mesh,ground,construction=fixture(4,False);base=construction['one_way']
    # Start from a valid MC reference ground, valid explicit mechanics and roles.
    jm=read(ROOT/'examples/joint-water-column-mesh.json');jc=read(ROOT/'examples/stages-joint-water-cycle.json')
    for policy in ('continuous','soil','prescribed'):
        b=copy.deepcopy(base);b['mechanics']['settings']['kinematics']='explicit';b['mechanics']['boundary']=jc['boundary'];b['mechanics']['joint_system']=copy.deepcopy(jc['joint_system'])
        config=b['mechanics']['joint_system'];config['joints'][0]['water']=dict(kind=policy)
        if policy=='soil':config['structure_bodies']=['UpperClay'];config['joints'][0]['water']['soil_side']='a'
        if policy=='prescribed':config['joints'][0]['water']['level_m']=.75
        # Only full-domain boundaries; JOINT rejection must precede any transfer.
        for period in b['flow']['analysis']['periods']:
            for bc in period['boundaries']:
                if bc['face_set']=='floor':bc['face_set']='top'
        check('ONE_WAY-'+policy,'solve-one-way',jm,ground,b,'ONE_WAY requires dry-reference weights')
        c=copy.deepcopy(construction);c['one_way']=copy.deepcopy(b)
        check('construction-JOINT-'+policy,'solve-one-way-construction',jm,ground,c,'ONE_WAY requires dry-reference weights')
    cases=[]
    b=copy.deepcopy(base);b['schema']='soil-fem-two-way/1';cases.append(('two-way-schema',b,'invalid ONE_WAY closure'))
    b=copy.deepcopy(base);b['coupling']='two_way';cases.append(('two-way-field',b,'unknown field'))
    b=copy.deepcopy(base);b['porosity_policy']='strain_feedback';cases.append(('porosity-feedback',b,'unknown variant'))
    b=copy.deepcopy(base);b['flow'].pop('initial_head');cases.append(('undefined-initial-head',b,'missing field'))
    b=copy.deepcopy(base);b['flow']['initial_head']={'kind':'previous_stage'};cases.append(('borrowed-initial-head',b,'unknown field'))
    b=copy.deepcopy(base);b['flow']['analysis']['periods'][1]['end_time_s']=.5;cases.append(('nonincreasing-physical-time',b,'period names/times'))
    b=copy.deepcopy(base);bc=b['flow']['analysis']['periods'][0]['boundaries'][0];bc['condition']['pressure_head_m']=.3;cases.append(('head-kind-value-mismatch',b,'unknown field'))
    b=copy.deepcopy(base);b['flow']['analysis']['periods'][0]['boundaries'].append(copy.deepcopy(b['flow']['analysis']['periods'][0]['boundaries'][0]));cases.append(('overlapping-heads',b,'overlapping/internal Richards boundary'))
    b=copy.deepcopy(base);b['flow']['joint_head_mapping']='incident_nodal_mean';cases.append(('fake-joint-head-map',b,'unknown field'))
    b=copy.deepcopy(base);b['flow']['time_s']=1.;cases.append(('fake-start-time',b,'unknown field'))
    # Use the soil-only controller on a full-domain valid period schedule.
    for name,b,expected in cases:
        for period in b['flow']['analysis']['periods']:
            for bc in period['boundaries']:
                if bc['face_set']=='floor':bc['face_set']='top'
        check(name,'solve-one-way',mesh,ground,b,expected)
    result=dict(passed=False,a04_complete=False,actual_cli_refusals=rows,
        bidirectional=dict(supported=False,controller=False,schema=False,cli=False,physical_solution_executed=False,
            evidence=['src/one_way.rs::Case and validate','src/one_way_construction.rs::Case','src/stages.rs::ExternalMechanics::new','src/main.rs command dispatch'],
            note='No coupling controller was added. Unknown schema/field refusals alone do not constitute a two-way solve.'),
        undefined_initial_head_not_filled=True,previous_stage_fields_not_borrowed=True)
    write(out/'refusals.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);a=p.parse_args();run(a.output.resolve(),a.solver.resolve())
