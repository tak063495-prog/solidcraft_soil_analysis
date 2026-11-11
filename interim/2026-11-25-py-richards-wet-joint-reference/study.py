"""Small joint-specific matrix; larger grids require explicit --pc."""
import argparse
import csv
import json
from pathlib import Path
import time
import numpy as np
from reference import Model,write,PARAMETERS,EVENTS
from compare import packet,NUMERIC


def difference(a,b):
    # Spatial comparison uses conservative restriction of fine CELL values, not
    # interpolation of joint traces. Different GP sets are never silently matched.
    ca=a['configuration'];cb=b['configuration']
    aa=np.asarray(a['final_head'])[:ca['nx']*ca['ny']].reshape(ca['ny'],ca['nx'])
    bb=np.asarray(b['final_head'])[:cb['nx']*cb['ny']].reshape(cb['ny'],cb['nx'])
    ry=cb['ny']//ca['ny'];rx=cb['nx']//ca['nx']
    if ry<1 or rx<1 or ca['ny']*ry!=cb['ny'] or ca['nx']*rx!=cb['nx']:raise ValueError('nested grid required')
    bb=bb.reshape(ca['ny'],ry,ca['nx'],rx).mean(axis=(1,3))
    return float(np.sqrt(np.mean((aa-bb)**2)))


def order(ds):
    return float(np.log2(ds[0]/ds[1])) if min(ds)>1e-14 else None


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--pc',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    report=dict(accepted=None,A04='未達',C04='未達',parameters=PARAMETERS,series={},cases={})
    def run(name,model,dt,end=1.):
        r=model.run(dt,end);write(a.output/(name+'.json'),r);return r
    try:
        for law in ('vg','linear'):
            for ny in (1,3):
                for mode,policy in [('ONE_WAY','CONTINUOUS'),('ONE_WAY','PRESCRIBED'),('TWO_WAY','CONTINUOUS')]:
                    name=f'{law}-{1 if ny==1 else 2}d-{mode}-{policy}'
                    r=run(name,Model(8,ny,law,mode,policy),.025)
                    report['cases'][name]={k:r[k] for k in ('maximum_water_balance_defect_m3','maximum_residual_m3_s','cumulative_net_inflow_m3')}
                    report['cases'][name].update(max_gap_m=max(x['gap_m'] for x in r['rows']),
                         max_joint_flow_m3_s=max(max(abs(x['q_a_to_joint_m3_s']),abs(x['q_joint_to_b_m3_s'])) for x in r['rows']),
                         contact_states=sorted(set(x['contact_state'] for x in r['rows'])))
                    if law=='vg' and ny==1 and mode=='ONE_WAY' and policy=='CONTINUOUS':
                        write(a.output/'reference-packet.json',packet(r))
        for mode in ('ONE_WAY','TWO_WAY'):
            r=run('dry-'+mode,Model(8,1,'vg',mode,dry=True),.025)
            report['cases']['dry-'+mode]=dict(maximum_water_balance_defect_m3=r['maximum_water_balance_defect_m3'],
              maximum_effective_pressure_kpa=max(x['effective_pressure_kpa'] for x in r['rows']),
              max_gap_m=max(x['gap_m'] for x in r['rows']))
        # Time and grid series are distinct. Self differences are not exact error.
        for law in ('vg','linear'):
            for mode in ('ONE_WAY','TWO_WAY'):
                temporal=[run(f'time-{law}-{mode}-{i}',Model(12,1,law,mode),dt) for i,dt in enumerate((.05,.025,.0125) if not a.pc else (.00625,.003125,.0015625))]
                ds=[difference(temporal[i],temporal[i+1]) for i in (0,1)]
                tables=[{(r['time_s'],r['point_id']):r for r in result['rows']} for result in temporal]
                common=sorted(set(tables[0]) & set(tables[1]) & set(tables[2]))
                interface={}
                for field in NUMERIC:
                    changes=[float(np.sqrt(np.mean([(tables[i][key][field]-tables[i+1][key][field])**2 for key in common]))) for i in (0,1)]
                    interface[field]=dict(rms_changes=changes,observed_order=order(changes))
                report['series'][f'time-{law}-{mode}']=dict(dt_s=[r['configuration']['dt_s'] for r in temporal],differences_rms_m=ds,observed_order=order(ds),interface_common_samples=len(common),interface=interface,contact_switch_brackets=[r['contact_switch_brackets'] for r in temporal])
            spatial=[run(f'space-{law}-{n}',Model(n,1,law,'ONE_WAY'),.003125 if not a.pc else .00078125) for n in ((4,8,16) if not a.pc else (32,64,128))]
            ds=[difference(spatial[i],spatial[i+1]) for i in (0,1)]
            report['series'][f'space-{law}']=dict(nx=[r['configuration']['nx'] for r in spatial],dt_s=spatial[0]['configuration']['dt_s'],differences_rms_m=ds,observed_order=order(ds),restriction='explicit pairwise cell-volume mean; not trace matching')
        # Opening-interface two-dimensional grid series (volume restriction only).
        for law in ('vg','linear'):
            grids=((4,2),(8,4),(16,8)) if not a.pc else ((16,8),(32,16),(64,32))
            results=[run(f'space2d-{law}-{nx}-{ny}',Model(nx,ny,law),.00625 if not a.pc else .0015625) for nx,ny in grids]
            ds=[difference(results[i],results[i+1]) for i in (0,1)]
            report['series'][f'space2d-{law}']=dict(grids=grids,dt_s=results[0]['configuration']['dt_s'],differences_rms_m=ds,observed_order=order(ds),restriction='explicit cell-volume means; no matching of different interface points')
        # Saturated exact eigenmode verifies the new interface assembly, independent
        # of C-04 RT0 and PR19 fixtures. Continuous, closed, zero gap storage.
        for nyfactor in (0,1):
            errors=[]
            for n in (4,8,16):
                ny=1 if nyfactor==0 else n//2
                dt=.0005 if nyfactor==0 else .00025
                r=run(f'eigen-{nyfactor}-{n}',Model(n,ny,'linear',schedule=False,eigen=True),dt,.02)
                z=np.repeat((np.arange(ny)+.5)/ny,n) if ny>1 else np.full(n,.25)
                xx=np.tile((np.arange(n)+.5)/n,ny)
                exact=2+.1*np.sin(np.pi*xx)*(np.cos(np.pi*z) if ny>1 else 1)*np.exp(-(2 if ny>1 else 1)*np.pi**2*.02)
                errors.append(float(np.sqrt(np.mean((np.array(r['final_head'][:n*ny])-exact)**2))))
            report['series'][f'eigen-{nyfactor}']=dict(errors_rms_m=errors,orders=[float(np.log2(errors[i]/errors[i+1])) for i in (0,1)],dt_s=dt,
              caveat='fixed time error remains; not pure asymptotic space order')
        base=json.loads((a.output/'vg-1d-ONE_WAY-CONTINUOUS.json').read_text())
        shifted=[]
        for shift in (-.0125,.0125):
            r=run(f'shift-{shift}',Model(8,1,switch_shift=shift),.025)
            shifted.append(dict(shift_s=shift,final_cell_head_difference_rms_m=difference(base,r),net_inflow_difference_m3=r['cumulative_net_inflow_m3']-base['cumulative_net_inflow_m3']))
        report['switch_sensitivity']=shifted
        # Shared exact samples only: no pressure/load-progress time alignment.
        ow=json.loads((a.output/'vg-1d-ONE_WAY-CONTINUOUS.json').read_text())
        tw=json.loads((a.output/'vg-1d-TWO_WAY-CONTINUOUS.json').read_text())
        report['coupling_difference']={k:float(np.sqrt(np.mean([(b[k]-aa[k])**2 for aa,b in zip(ow['rows'],tw['rows'],strict=True)]))) for k in NUMERIC}
        report.update(complete=True,elapsed_seconds=time.monotonic()-started)
        write(a.output/'summary.json',report)
    except Exception as error:
        write(a.output/'failure.json',dict(complete=False,accepted=None,error=repr(error)));raise

if __name__=='__main__':main()
