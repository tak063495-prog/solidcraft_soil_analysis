"""Isolated curved-geometry guards, typed declarations and legacy CLI controls."""
import argparse,copy,sys
from pathlib import Path
from b03_cases import ROOT,make,detached_fold,write,read
from study import command,sha


def run(out,cli):
    out.mkdir(parents=True,exist_ok=False);rows=[]
    def record(name,data,expected=None,legacy=False):
        m,g,c,_=data;f=out/name;f.mkdir()
        for k,v in [('mesh',m),('ground',g),('stages',c)]:write(f/(k+'.json'),v)
        target=f/'run.json';p=command([cli,'solve-stages',*(f/(k+'.json') for k in ('mesh','ground','stages')),target],f/'cli.log',expect=None)
        log=(f/'cli.log').read_text()
        row=dict(name=name,expected=expected,legacy_control=legacy,process=p,diagnostic=log,result_exists=target.exists(),
                 complete=read(target).get('complete') if target.exists() else None,
                 refused=p['exit_code']==1 and expected is not None and expected in log and not target.exists(),
                 input_sha256={k:sha(f/(k+'.json')) for k in ('mesh','ground','stages')})
        rows.append(row);write(out/'supplemental.json',dict(passed=False,b03_complete=False,records=rows))
    for axis in (0,2):
        data=detached_fold(0,0,0);cfg=data[2]['beam_system'];n=next(n for n in cfg['nodes'] if n['id']==cfg['plates'][0]['nodes'][3]);n['position_m'][axis]+=.001
        record('isolated-curved-midpoint-'+str(axis),data,'Mindlin Tri6 requires straight edges')
    data=detached_fold(90,.05,.05)
    record('implicit-planar-formulation-with-real-fold',data,'same signed normal offset')
    for kind in ('mindlin_tri6','cst_dkt','default'):
        data=make('kink-zero');cfg=data[2]['beam_system']
        for p in cfg['plates']:
            p.pop('attachment_offset_m')
            if kind!='default':p['formulation']=kind
            else:p.pop('formulation')
            if kind!='mindlin_tri6':p['nodes']=p['nodes'][:3]
        record('legacy-'+kind+'-kink-zero',data,legacy=True)
    # Execute the existing compact entry point as a separate Python process.
    m,g,c,_=make('planar');base=out/'compact';base.mkdir();c.pop('beam_system')
    spec=dict(schema='soil-fem-plate-surfaces/1',sections=[dict(id='wall',young_kpa=3e7,poisson=.2,thickness_m=.1,unit_weight_kn_m3=24.,drilling_factor=1e-4)],
              surfaces=[dict(id='top',section_id='wall',face_set='top',initially_active=False,formulation='mindlin_tri6',thickness_offset_m=.05)])
    for formulation,offset,expected in [('mindlin_tri6',.05,'plate thickness offset attachment unsupported'),('mindlin_tri6_offset',0.,'unknown plate formulation')]:
        name=f'{formulation}-{offset}';f=base/name;f.mkdir();spec['surfaces'][0].update(formulation=formulation,thickness_offset_m=offset)
        for k,v in [('mesh',m),('stages',c),('specification',spec)]:write(f/(k+'.json'),v)
        p=command([sys.executable,ROOT/'tools/plate_surfaces.py',f/'mesh.json',f/'stages.json',f/'specification.json',f/'expanded'],f/'cli.log',expect=None)
        log=(f/'cli.log').read_text();rows.append(dict(name='compact-'+name,process=p,diagnostic=log,expected=expected,refused=p['exit_code']!=0 and expected in log and not (f/'expanded').exists()))
    # Native preparation guard reached directly; no GUI or CAD claim.
    sys.path.insert(0,str(ROOT/'tools'));from prepare_native_plates import derive
    for kind in ('planar','kink-zero'):
        m,g,c,_=make(kind)
        try:derive(m,c)
        except ValueError as e:rows.append(dict(name='native-prepare-'+kind,entry='prepare_native_plates.derive',diagnostic=str(e),refused='Native offset plate' in str(e)))
        else:rows.append(dict(name='native-prepare-'+kind,refused=False))
    write(out/'supplemental.json',dict(passed=False,b03_complete=False,records=rows,
         refusal_contract_passed=all(r.get('refused') for r in rows if not r.get('legacy_control')),
         legacy_controls_completed=all(r['process']['exit_code']==0 and r['complete'] for r in rows if r.get('legacy_control'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--cli',type=Path,required=True);a=p.parse_args();run(a.output.resolve(),a.cli.resolve())
