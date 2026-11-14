"""Actual production CLI accept/refuse matrix, never substitute geometry audit for solve."""
import argparse,copy
from pathlib import Path
import numpy as np
from b03_cases import make,detached_fold,write,B02
from study import command,sha


def run(out,cli):
    out.mkdir(parents=True,exist_ok=False);rows=[]
    def check(name,data,expected=None,category='geometry'):
        m,g,c,_=data;f=out/name;f.mkdir()
        for k,v in [('mesh',m),('ground',g),('stages',c)]:write(f/(k+'.json'),v)
        target=f/'run.json';meta=command([cli,'solve-stages',*(f/(k+'.json') for k in ('mesh','ground','stages')),target],f/'cli.log',expect=None)
        diagnostic=(f/'cli.log').read_text();exists=target.exists()
        row=dict(name=name,category=category,expected_diagnostic=expected,diagnostic=diagnostic,process=meta,
                 result_file_exists=exists,complete=__import__('json').loads(target.read_text()).get('complete') if exists else None,
                 input_sha256={k:sha(f/(k+'.json')) for k in ('mesh','ground','stages')})
        row['refused_as_expected']=meta['exit_code']==1 and not exists and expected is not None and expected in diagnostic
        row['accepted']=meta['exit_code']==0 and row['complete'] is True
        rows.append(row);write(out/'scope.json',dict(passed=False,b03_complete=False,records=rows))
        print(name,'accepted' if row['accepted'] else diagnostic.strip()[:160],flush=True)
    for kind in ('planar','kink-zero','normal-sign-pair','step-split','kink-split','height-step'):
        check(kind,make(kind),category='accepted_control')
    for angle in (15,45,90,135,179):
        for tag,e1,e2 in [('positive',.05,.05),('negative',-.05,-.05),('one-zero',.05,0.)]:
            check(f'real-fold-{angle}-{tag}',detached_fold(angle,e1,e2),'same signed normal offset')
    def variant(name,fn,kind='planar',expected='same signed normal offset',category='mutation'):
        data=copy.deepcopy(make(kind));fn(*data[:3]);check(name,data,expected,category)
    def plates(c):return c['beam_system']['plates']
    variant('step-e',lambda m,g,c:plates(c)[1].update(attachment_offset_m=.04))
    variant('sign-only',lambda m,g,c:plates(c)[1].update(attachment_offset_m=-.05))
    def reverse(m,g,c):
        p=plates(c)[1];p['nodes']=[p['nodes'][j] for j in (0,2,1,5,4,3)]
    variant('normal-only',reverse)
    for label,e in [('positive',.051),('negative',-.051)]:
        variant('outside-'+label,lambda m,g,c,e=e:plates(c)[0].update(attachment_offset_m=e),expected='finite |e| <= thickness/2')
    def curve(m,g,c,axis=2,structure_only=False):
        name=plates(c)[0]['nodes'][3];n=next(n for n in c['beam_system']['nodes'] if n['id']==name)
        n['position_m']=n['position_m'].copy();n['position_m'][axis]+=.001
        if not structure_only:m['nodes'][m['node_sets'][n['soil_node_set']][0]]=n['position_m'].copy()
    variant('curved-midpoint-out-of-plane',curve,expected='Mindlin Tri6')
    variant('curved-midpoint-in-plane',lambda m,g,c:curve(m,g,c,0),expected='Mindlin Tri6')
    variant('one-sided-attachment-change',lambda m,g,c:curve(m,g,c,structure_only=True),expected='attachment coordinate mismatch')
    variant('claimed-planar-real-fold',lambda m,g,c:c['beam_system'].update(reference_surface='planar'),kind='kink-zero',expected='unknown field',category='typed')
    variant('per-node-normal',lambda m,g,c:c['beam_system']['nodes'][0].update(normal=[0,0,1]),expected='unknown field',category='typed')
    variant('per-node-e',lambda m,g,c:c['beam_system']['nodes'][0].update(attachment_offset_m=.05),expected='unknown field',category='typed')
    variant('old-mindlin-nonzero-e',lambda m,g,c:plates(c)[0].update(formulation='mindlin_tri6'),expected='requires explicit mindlin_tri6_offset')
    variant('old-cst-nonzero-e',lambda m,g,c:plates(c)[0].update(formulation='cst_dkt',nodes=plates(c)[0]['nodes'][:3]),expected='requires explicit mindlin_tri6_offset')
    result=dict(passed=False,b03_complete=False,records=rows,
                refusal_contract_passed=all(r['refused_as_expected'] for r in rows if r['expected_diagnostic'] is not None),
                accepted_controls_passed=all(r['accepted'] for r in rows if r['category']=='accepted_control'))
    write(out/'scope.json',result)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--cli',type=Path,required=True);a=p.parse_args();run(a.output.resolve(),a.cli.resolve())
