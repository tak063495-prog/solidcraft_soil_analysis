"""Single-body planar-surface example seed; equilibrium is NOT assumed.

The input surface elevation, K0, and pc/p are explicit. No pressure clipping is
used. The generated per-point stresses must pass solver admissibility checks
and geostatic iteration. General imported states need not use this generator.
"""
import argparse
import hashlib
import json
from pathlib import Path
import math


def create(mesh,identity,template,output,top,gx,gy,kx,ky,ratio):
    if output.exists():raise ValueError('seed input exists')
    m=json.loads(mesh.read_text());layout=json.loads(identity.read_text());g=json.loads(template.read_text())
    assert layout['schema']=='soil-fem-mesh-identity/1' and g['initialization']['method']=='equilibrate'
    if len(g['bodies'])!=1 or len(g['materials'])!=1:raise ValueError('this planar overburden generator supports a single body/material only')
    if not all(math.isfinite(v) for v in (top,gx,gy,kx,ky,ratio)) or min(kx,ky,ratio)<=0:raise ValueError('invalid seed parameters')
    material=g['materials'][0];cfg=g['initialization']['geostatic'];points=[]
    for p in layout['integration_points']:
        x,y,z=p['position_m'];depth=top+gx*x+gy*y-z
        total=g['initialization']['surface_surcharge_kpa']+material['unit_weight_kn_m3']*depth
        water=g['initialization']['water_level_m'];pore=0 if water is None else max(water-z,0)*g['initialization']['water_unit_weight_kn_m3']
        vertical=total-pore
        if depth<=0 or vertical<=0:raise ValueError('seed point has nonpositive depth/effective overburden; no clipping is applied')
        stress=[kx*vertical,ky*vertical,vertical,0.,0.,0.];state={'effective_stress_compression_positive_kpa':stress,'plastic_strain_compression_positive':[0.]*6}
        if 'mcc' in material:
            a=material['mcc'];mean=sum(stress[:3])/3;pc=ratio*mean
            q2=1.5*sum((s-mean)**2 for s in stress[:3]);yield_value=(q2+a['critical_state_m']**2*mean*(mean-pc))/(a['critical_state_m']**2*pc**2)
            if yield_value>1e-8:raise ValueError('requested K0/pc ratio puts the seed outside MCC')
            e=a['specific_volume_ncl_at_reference']-a['lambda']*math.log(pc/a['reference_pressure_kpa'])+a['kappa']*math.log(pc/mean)-1
            if e<=0:raise ValueError('MCC seed void ratio is nonpositive')
            state.update(preconsolidation_pressure_kpa=pc,void_ratio=e)
        else:
            seed=next(iter(cfg['seed']['states'].values()));state.update(void_ratio=seed['void_ratio'],equivalent_plastic_strain=0.,plastic_work_kpa=0.,last_return='elastic')
        points.append({'element_index':p['element_index'],'local_index':p['local_index'],'position_m':p['position_m'],'material_id':m['material_ids'][p['element_index']],'state':state})
    cfg['seed']={'kind':'integration_points','mesh_sha256':layout['mesh_sha256'],'points':points}
    output.write_text(json.dumps(g,allow_nan=False))
    record={'schema':'soil-fem-planar-overburden-seed/1','equilibrium_assumed':False,'surface_z':'top + gradient_x*x + gradient_y*y',
        'top_at_origin_m':top,'gradient_x':gx,'gradient_y':gy,'k0_x':kx,'k0_y':ky,'mcc_pc_over_mean_effective':ratio,
        'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (mesh,identity,template,output)}}
    output.with_suffix('.seed-provenance.json').write_text(json.dumps(record,indent=2));print(f'Generated {len(points)} unbalanced initial stress points: {output}')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('identity',type=Path);p.add_argument('template',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--top-at-origin-m',type=float,required=True);p.add_argument('--gradient-x',type=float,default=0);p.add_argument('--gradient-y',type=float,default=0)
    p.add_argument('--k0-x',type=float,required=True);p.add_argument('--k0-y',type=float,required=True);p.add_argument('--pc-ratio',type=float,required=True)
    a=p.parse_args();create(a.mesh,a.identity,a.template,a.output,a.top_at_origin_m,a.gradient_x,a.gradient_y,a.k0_x,a.k0_y,a.pc_ratio)
