"""Streaming common-physical-coordinate differences; no extrapolated fields."""
import argparse
import math
from pathlib import Path
import numpy as np
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import mapping as m
c=m.c
FIELDS=m.NODAL+m.MATERIAL+m.CONTACT+('open','void_from_strain','pore_physical')

def order_metric(a,b,dot):
    return dict(coarse_absolute_l2=float(np.sqrt(a)),fine_absolute_l2=float(np.sqrt(b)),observed_order=float(np.log2(np.sqrt(a/b))) if a*b>0 else None,difference_direction_cosine=float(dot/np.sqrt(a*b)) if a*b>0 else None)

def selected(out,h,g,kind,n):
    assert g%n==0
    count=0
    for i,(row,f) in enumerate(c.old.states(out/f'h{h}-g{g}'/kind)):
        if (i%g+1)%(g//n)==0:
            phase,step=divmod(count,n)
            assert row['stage']==c.STAGES[phase] and row['progress']==(step+1)/n
            yield row,f;count+=1
    assert count==7*n

class Stats:
    def __init__(self):self.values={};self.stage={}
    def add(self,a,b,stage):
        for k in a:
            d=(a[k]-b[k]).reshape(len(a[k]),-1);r=b[k].reshape(len(b[k]),-1)
            d2=(d*d).sum(axis=1);r2=(r*r).sum(axis=1);mx=abs(d).max(axis=1)
            if k not in self.values:self.values[k]=[np.zeros_like(d2),np.zeros_like(r2),np.zeros_like(mx)]
            v=self.values[k];v[0]+=d2;v[1]+=r2;v[2]=np.maximum(v[2],mx)
            s=self.stage.setdefault(stage,{}).setdefault(k,[0.,0.,0.]);s[0]+=float(d2.sum());s[1]+=float(r2.sum());s[2]=max(s[2],float(mx.max()))
    @staticmethod
    def metric(d2,r2,mx):
        return dict(absolute_difference_l2=float(np.sqrt(d2)),difference_squared=float(d2),reference_squared=float(r2),pooled_relative_l2=float(np.sqrt(d2)/max(np.sqrt(r2),1e-12)),maximum_absolute=float(mx))
    def report(self):
        return dict(fields={k:self.metric(v[0].sum(),v[1].sum(),v[2].max()) for k,v in self.values.items()},by_stage={s:{k:self.metric(*v) for k,v in fs.items()} for s,fs in self.stage.items()})
    def save(self,path):
        with path.open('xb') as f:np.savez_compressed(f,**{k+'__'+suffix:v[i] for k,v in self.values.items() for i,suffix in enumerate(('difference_squared','reference_squared','maximum_absolute'))})
        return dict(file=path.name,sha256=c.sha(path),entity_counts={k:len(v[0]) for k,v in self.values.items()})

def triple(out,target,maps,routes,kind,label):
    n=math.gcd(*(g for h,g in routes));iterators=[selected(out,h,g,kind,n) for h,g in routes]
    stats=[Stats(),Stats()];order={};sensitivity=[Stats() for _ in routes];defects=[{} for _ in routes]
    status_stats=[Stats(),Stats()];count=0
    enums=[]
    for h,g in routes:
        p=out/f'h{h}-g{g}'/'endpoint-states.json';enums.append(c.read(p) if kind=='native' else None)
    for states in zip(*iterators,strict=True):
        assert len({(r['stage'],r['progress']) for r,f in states})==1
        stage=states[0][0]['stage'];arrays=[];statuses=[]
        for j,((h,g),(row,raw)) in enumerate(zip(routes,states)):
            a=maps.apply(h,raw);arrays.append(a)
            sensitivity[j].add(a,maps.apply(h,raw,True),stage)
            # Known analytic branch: interpolation error of exponential void and
            # cut pressure can be measured independently of the h difference.
            strain=a['strain'];owners=maps.operators[h]['material_owners'];mesh=maps.meshes[h]
            region={r['id']:r['body'] for r in mesh['regions']};seeds=c.read(out/'ground.json')['initialization']['geostatic']['seed']['states']
            e0=np.array([seeds[region[mesh['region_ids'][e]]]['void_ratio'] for e in owners])
            expected_void=(1+e0)*np.exp(strain[:,:3].sum(axis=1))-1
            ground=c.read(out/'ground.json');expected_pore=ground['initialization']['water_unit_weight_kn_m3']*np.maximum(ground['initialization']['water_level_m']-maps.material_position[:,2],0)
            for k,d in [('void_analytic_interpolation_error',a['void']-expected_void),('pore_affine_cut_interpolation_error',a['pore']-expected_pore)]:
                defects[j][k]=max(defects[j].get(k,0.),float(abs(d).max()))
            normal=maps.systems[h].normal[0]
            assert np.max(abs(maps.systems[h].normal-normal))<1e-12
            qn=a['traction']@normal
            kn=maps.systems[h].kn[0];mu=maps.systems[h].mu[0]
            tau=a['traction']-qn[:,None]*normal
            extra=[('contact_normal_law_reconstruction_defect_kpa',qn-kn*np.minimum(a['gap'],0)),
                   ('contact_open_traction_reconstruction_defect_kpa',np.linalg.norm(a['traction'][a['gap']>0],axis=1)),
                   ('contact_friction_reconstruction_excess_kpa',np.maximum(np.linalg.norm(tau,axis=1)-mu*np.maximum(-qn,0),0))]
            for key,d in extra:defects[j][key]=max(defects[j].get(key,0.),float(np.max(abs(d),initial=0)))
            if enums[j] is not None:
                raw_index=c.STAGES.index(stage)*g+round(row['progress']*g)-1
                codes={'open':0.,'stick':1.,'slip':2.}
                status=np.array([codes[s] for s in enums[j][raw_index]])
                # Categorical nearest point, never interpolate enum integers.
                statuses.append(maps.operators[h]['nearest_contact']@status)
        for j in range(2):
            stats[j].add(arrays[j],arrays[j+1],stage)
            if statuses:
                mismatch=(statuses[j]!=statuses[j+1]).astype(float)
                status_stats[j].add({'status_mismatch':mismatch},{'status_mismatch':np.zeros_like(mismatch)},stage)
        for k in FIELDS:
            d0=arrays[0][k]-arrays[1][k];d1=arrays[1][k]-arrays[2][k]
            v=order.setdefault(k,[0.,0.,0.]);v[0]+=float((d0*d0).sum());v[1]+=float((d1*d1).sum());v[2]+=float((d0*d1).sum())
        count+=1
    result=dict(routes=routes,kind=kind,common_states=count,time_interpolation=False,all_physical_coordinates_identical=True,pairs=[],orders={})
    for j,s in enumerate(stats):
        report=s.report();report['per_entity']=s.save(target/f'{label}-pair{j}-entities.npz')
        if statuses:report['categorical_nearest_status']=status_stats[j].report()
        result['pairs'].append(report)
    for k,(a,b,dot) in order.items():
        result['orders'][k]=order_metric(a,b,dot)
    result['reconstruction_sensitivity']=[s.report() for s in sensitivity];result['known_branch_interpolation_defects']=defects
    result.update(accepted=False,order_definition='log2(||x0-x1||/||x1-x2||), same fixed coordinates and weights; simultaneous series is a combined h/path slope')
    print('measured',label,flush=True);return result

def pair(out,target,maps,routes,kind,label):
    n=math.gcd(*(g for h,g in routes));stat=Stats();count=0
    for states in zip(*(selected(out,h,g,kind,n) for h,g in routes),strict=True):
        assert len({(r['stage'],r['progress']) for r,f in states})==1
        a,b=[maps.apply(h,f) for (h,g),(r,f) in zip(routes,states)]
        stat.add(a,b,states[0][0]['stage']);count+=1
    result=stat.report();result.update(routes=routes,kind=kind,common_states=count,per_entity=stat.save(target/f'{label}-entities.npz'),accepted=False)
    return result

def contamination(out,maps,h0,h1,ga,gb):
    routes=[(h0,ga),(h1,ga),(h0,gb),(h1,gb)];n=math.gcd(ga,gb);values={};count=0
    for states in zip(*(selected(out,h,g,'native',n) for h,g in routes),strict=True):
        arrays=[maps.apply(h,f) for (h,g),(row,f) in zip(routes,states)]
        for k in FIELDS:
            sa=arrays[0][k]-arrays[1][k];sb=arrays[2][k]-arrays[3][k]
            delta=sb-sa;ta=arrays[2][k]-arrays[0][k];tb=arrays[3][k]-arrays[1][k]
            assert np.allclose(delta,ta-tb,rtol=0,atol=1e-12)
            v=values.setdefault(k,[0.,0.,0.,0.,0.]);v[0]+=float((sa*sa).sum());v[1]+=float((sb*sb).sum());v[2]+=float((delta*delta).sum());v[3]+=float((sa*sb).sum());v[4]=max(v[4],float(abs(delta-(ta-tb)).max()))
        count+=1
    return dict(routes=routes,common_states=count,identity='S(gb)-S(ga) = T(h0)-T(h1)',fields={k:dict(spatial_ga_absolute_l2=float(np.sqrt(a)),spatial_gb_absolute_l2=float(np.sqrt(b)),change_of_spatial_difference_l2=float(np.sqrt(d)),change_over_gb_spatial_difference=float(np.sqrt(d/b)) if b>0 else None,spatial_difference_direction_cosine=float(dot/np.sqrt(a*b)) if a*b>0 else None,decomposition_identity_maximum_roundoff=err) for k,(a,b,d,dot,err) in values.items()},sufficiently_fine_time_not_certified=True)

def run(out,target,max_level=2):
    target.mkdir(parents=True,exist_ok=False)
    maps=m.Maps([c.read(out/f'mesh-h{h}.json') for h in range(max_level+1)],c.read(out/'stages.json'))
    c.write(target/'common-coordinates.json',maps.coordinates())
    report=dict(schema='wet-joint-common-physical-coordinate-convergence/1',starting_commit=c.BASE,passed=False,c02_judgement='未達',p6_complete=False,c01_c03_c04_accepted=False,spatial={},temporal_contamination={},implementation_comparison={},restarts={})
    series=[(f'fixed{g}',[(0,g),(1,g),(2,g)]) for g in (32,64,128,256,512)]+[('simultaneous',[(0,16),(1,32),(2,64)])]
    for label,routes in series:
        if all((out/f'h{h}-g{g}/native/manifest.json').exists() for h,g in routes):report['spatial'][label]=triple(out,target,maps,routes,'native',label)
    for h in range(max_level+1):
        for a,b in [(32,64),(64,128),(128,256),(256,512)]:
            routes=[(h,a),(h,b)]
            if all((out/f'h{h}-g{g}/native/manifest.json').exists() for h,g in routes):report['temporal_contamination'][f'h{h}-g{a}-{b}']=pair(out,target,maps,routes,'native',f'time-h{h}-{a}-{b}')
    for g in (32,64,128,256,512):
        routes=[(0,g),(1,g)]
        if all((out/f'h{h}-g{g}/native/manifest.json').exists() for h,g in routes):report['temporal_contamination'][f'spatial-h01-g{g}']=pair(out,target,maps,routes,'native',f'space-h01-{g}')
    report['contamination_decomposition']={}
    for h0,h1,ga,gb in [(0,1,32,64),(1,2,32,64),(0,1,64,256),(0,1,128,256)]:
        if all((out/f'h{h}-g{g}/native/manifest.json').exists() for h in (h0,h1) for g in (ga,gb)):
            report['contamination_decomposition'][f'h{h0}{h1}-g{ga}-{gb}']=contamination(out,maps,h0,h1,ga,gb)
    # Same-mesh comparison preserves every original GP/node; no reconstruction.
    from importlib.util import spec_from_file_location,module_from_spec
    spec=spec_from_file_location('c02_old_measure',c.ROOT/'interim/2026-11-15-wet-joint-convergence/measure.py');oldmeasure=module_from_spec(spec);spec.loader.exec_module(oldmeasure)
    for h in range(max_level+1):
        for g in (16,32,64,128,256,512):
            folder=out/f'h{h}-g{g}'
            if all((folder/k/'manifest.json').exists() for k in ('native','independent')):
                report['implementation_comparison'][f'h{h}-g{g}']=oldmeasure.compare(folder/'native',folder/'independent')
                oldmeasure.load.cache_clear()
            if (folder/'restart-verification.json').exists():report['restarts'][f'h{h}-g{g}']=c.read(folder/'restart-verification.json')
    for h in range(max_level-1):
        if h==0:continue
        for g in (32,64,128,256,512):
            routes=[(j,g) for j in (h,h+1,h+2)]
            if all((out/f'h{j}-g{g}/native/manifest.json').exists() for j,g in routes):
                label=f'fixed{g}-h{h}{h+1}{h+2}';report['spatial'][label]=triple(out,target,maps,routes,'native',label)
        routes=[(h,32),(h+1,64),(h+2,128)]
        if all((out/f'h{j}-g{g}/native/manifest.json').exists() for j,g in routes):
            label=f'simultaneous-h{h}{h+1}{h+2}';report['spatial'][label]=triple(out,target,maps,routes,'native',label)
    report['unmeasured_cloud_fine_fixed_time_routes']=[f'h2-g{g}' for g in (128,256) if not (out/f'h2-g{g}/native/manifest.json').exists()]
    report['limitations']=['Contact constitutive histories have no interpolation error bound between GPs; nearest/cubic sensitivity is not a bound.','Fine fixed steps must not be called sufficiently fine solely because N=256; C-01 remains unresolved.','No continuum reference, extrapolated field, correction, new tolerance, or general plastic/water acceptance.']
    c.write(target/'measurements.json',report)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--max-level',type=int,default=2);a=p.parse_args();run(a.input.resolve(),a.output.resolve(),a.max_level)
