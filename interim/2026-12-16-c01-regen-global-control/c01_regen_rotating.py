"""Manufactured two-tangent paths: input-only exact continuous solutions.

Known traction c(t)n(t), plastic rate gamma*n(t), and x'=tau'/kt+s'
construct each sliding arc. Elastic segments stay inside the Coulomb disk.
Opening unloads c to zero before an open excursion and fresh recontact.
No native output is used to construct these inputs or their exact integrals.
"""
import argparse
import json
from math import cos, sin, pi, fsum, log2, sqrt
from pathlib import Path
import subprocess
import numpy as np

KN, KT, MU = 1e6, 2e5, .4


def norm(v):
    return sqrt(fsum(float(x)*float(x) for x in np.asarray(v).flat))


def paths():
    # (kind, end cap, angle change, plastic distance), unit duration per segment.
    return {
        'arc_variable_pressure': [('elastic',40.,0.,0.),('arc',120.,pi, .004),('arc',60.,pi/2,.002)],
        'arc_reverse': [('elastic',80.,0.,0.),('arc',40.,pi,.003),('turn',40.,pi,0.),('arc',100.,-pi,.003)],
        'polygon_open_recontact': [('elastic',80.,0.,0.),('arc',40.,0.,.002),('turn',40.,pi/2,0.),('arc',100.,0.,.002),('unload',0.,0.,0.),('open',0.,0.,0.),('elastic',60.,0.,0.),('arc',30.,-pi/2,.003)],
    }


def build(name, outer):
    segments=paths()[name]
    jump=np.zeros(3);plastic=np.zeros(3);tau=np.zeros(2);d=0.;theta=0.;cap=0.
    samples=[]
    for si,(kind,endcap,angle,gamma) in enumerate(segments):
        x0=jump.copy();p0=plastic.copy();q0=tau.copy();d0=d;c0=cap;theta0=theta
        for k in range(1,outer+1):
            t=k/outer;c=c0+(endcap-c0)*t
            if kind=='open':
                # A triangular positive gap, with nonzero tangential excursion;
                # all grids include its midpoint; no traction or plastic flow.
                q=np.zeros(2);dp=np.zeros(2)
                x=x0[:2]+np.array([.001*t,-.0007*t]);g=.0002*(2*t if t<=.5 else 2*(1-t))
                exact_d=d0
            elif kind in ('elastic','unload'):
                direction=np.array([cos(theta0),sin(theta0)])
                q=c*direction;dp=np.zeros(2);x=x0[:2]+(q-q0)/KT;g=-c/(MU*KN);exact_d=d0
            elif kind=='turn':
                end=endcap*np.array([cos(theta0+angle),sin(theta0+angle)])
                q=q0+(end-q0)*t;dp=np.zeros(2);x=x0[:2]+(q-q0)/KT;g=-c/(MU*KN);exact_d=d0
            else:
                a=theta0+angle*t;q=c*np.array([cos(a),sin(a)])
                dp=(gamma*t*np.array([cos(theta0),sin(theta0)]) if angle==0 else gamma/angle*np.array([sin(a)-sin(theta0),cos(theta0)-cos(a)]))
                x=x0[:2]+(q-q0)/KT+dp;g=-c/(MU*KN)
                exact_d=d0+gamma*(c0*t+(endcap-c0)*t*t/2)
            dest=np.r_[x,g];ps=p0+np.r_[dp,0.]
            samples.append(dict(segment=si,coordinate=t,jump=dest.tolist(),plastic_slip=ps.tolist(),traction=[*q.tolist(),KN*min(g,0.)],dissipation=exact_d,sliding=kind=='arc',direction=[cos(theta0+angle*t),sin(theta0+angle*t),0.]))
        jump=dest;plastic=ps;tau=q;d=exact_d;cap=endcap;theta+=angle
    return samples


def request(outers=(2,4,8,16,32,64),counts=(1,2,4,8)):
    runs=[];truth={}
    for name in paths():
        for n in outers:
            samples=build(name,n)
            for count in counts:
                key=f'{name}/n{n}/k{count}'
                runs.append(dict(id=key,policy='fixed',fixed_count=count,steps=[dict(jump_m=r['jump'],trial_jumps_m=[]) for r in samples]))
                truth[key]=samples
    return dict(material=dict(normal_stiffness_kpa_m=KN,tangential_stiffness_kpa_m=KT,friction_coefficient=MU),normal=[0.,0.,1.],maximum_increment_m=5e-5,runs=runs),truth


def compare(response,truth):
    # Independent discrete vector update; independent analytic manufactured
    # continuous solution above is the precision reference, with zero truncation.
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'2026-12-15-c01-local-path-decomposition'))
    from c01_path_replay import integrate
    rows={}
    for run in response:
        p=np.zeros((1,3));q=p.copy();x=p.copy();d=0.;old_native_d=old_exact_d=0.
        errors={k:[] for k in ('slip','traction','increment_D','cumulative_D','endpoint_implementation','direction_radians','slip_increment_direction_radians')};history=[];old_ns=np.zeros(3);old_ps=np.zeros(3)
        for sample,exact in zip(run['samples'],truth[run['id']],strict=True):
            dest=np.array([exact['jump']]);qq,dp,dd,_,_=integrate(x,q,dest,np.array([[0.,0.,1.]]),np.array([KN]),np.array([KT]),np.array([MU]),sample['actual_count'])
            p+=dp;d+=float(dd[0]);q=qq;x=dest
            state=sample['state'];ns=np.array(state['plastic_slip_m']);nt=np.array(state['traction_kpa']);nd=state['plastic_dissipation_kpa_m']
            errors['slip'].append(float(norm(ns-exact['plastic_slip'])))
            errors['traction'].append(float(norm(nt-exact['traction'])))
            errors['increment_D'].append(abs((nd-old_native_d)-(exact['dissipation']-old_exact_d)))
            errors['cumulative_D'].append(abs(nd-exact['dissipation']))
            errors['endpoint_implementation'].append(max(float(np.max(abs(ns-p[0]))),float(np.max(abs(nt-q[0]))),abs(nd-d)))
            if exact['sliding'] and norm(nt[:2])>0:
                # Angle against instantaneous exact plastic direction, not a
                # previous step direction or total chord displacement.
                a=np.arctan2(nt[1],nt[0]);b=np.arctan2(exact['direction'][1],exact['direction'][0]);errors['direction_radians'].append(abs(float(np.arctan2(np.sin(a-b),np.cos(a-b)))))
            nsinc=ns-old_ns;psinc=np.asarray(exact['plastic_slip'])-old_ps
            if exact['sliding'] and norm(nsinc)>0 and norm(psinc)>0:
                a=np.arctan2(nsinc[1],nsinc[0]);b=np.arctan2(psinc[1],psinc[0]);errors['slip_increment_direction_radians'].append(abs(float(np.arctan2(np.sin(a-b),np.cos(a-b)))))
            old_ns=ns;old_ps=np.asarray(exact['plastic_slip'])
            history.append(dict(exact=exact,native=sample,independent_discrete_slip=p[0].tolist(),independent_discrete_traction=q[0].tolist(),independent_discrete_D=d))
            old_native_d=nd;old_exact_d=exact['dissipation']
        rows[run['id']]=dict(maximum_errors={k:max(v,default=0.) for k,v in errors.items()},minimum_increment_D=min(s['plastic_dissipation_increment_kpa_m'] for s in run['samples']),final_D_error=nd-truth[run['id']][-1]['dissipation'],history=history)
    series={}
    for name in paths():
        for count in (1,2,4,8):
            keys=[f'{name}/n{n}/k{count}' for n in (2,4,8,16,32,64)]
            if not all(k in rows for k in keys):continue
            fields={}
            for field in ('slip','traction','increment_D','cumulative_D','direction_radians','slip_increment_direction_radians'):
                es=[rows[k]['maximum_errors'][field] for k in keys]
                fields[field]=dict(errors=es,monotone=all(b<=a for a,b in zip(es,es[1:])),orders=[log2(a/b) if a*b>1e-28 else None for a,b in zip(es,es[1:])])
            series[f'{name}/k{count}']=fields
    return dict(schema='c01-regen-rotating/1',c01_passed=False,exact_reference='manufactured Coulomb spring solution; analytic rotating plastic-rate integral, linear cap integral; elastic disk chords; zero-traction open excursion',runs=rows,series=series)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    req,truth=request();(a.output/'request.json').write_text(json.dumps(req)+'\n')
    subprocess.run([str(a.solver.resolve()),str((a.output/'request.json').resolve()),str((a.output/'response.json').resolve())],check=True)
    response=json.loads((a.output/'response.json').read_text())
    result=compare(response,truth)
    from copy import deepcopy
    rz=np.array([[cos(.37),-sin(.37),0.],[sin(.37),cos(.37),0.],[0.,0.,1.]])
    rx=np.array([[1.,0.,0.],[0.,cos(.71),-sin(.71)],[0.,sin(.71),cos(.71)]])
    rotation=rx@rz;rotated=deepcopy(req);rotated['normal']=(rotation@np.array(req['normal'])).tolist()
    for run in rotated['runs']:
        for step in run['steps']:step['jump_m']=(rotation@np.array(step['jump_m'])).tolist()
    (a.output/'rotation-request.json').write_text(json.dumps(rotated)+'\n')
    subprocess.run([str(a.solver.resolve()),str((a.output/'rotation-request.json').resolve()),str((a.output/'rotation-response.json').resolve())],check=True)
    transformed=json.loads((a.output/'rotation-response.json').read_text());errors={k:[] for k in ('plastic_slip_m','traction_kpa','plastic_dissipation_kpa_m')}
    for ordinary,changed in zip(response,transformed,strict=True):
        for left,right in zip(ordinary['samples'],changed['samples'],strict=True):
            for key in ('plastic_slip_m','traction_kpa'):
                errors[key].append(norm(np.array(right['state'][key])-rotation@np.array(left['state'][key])))
            errors['plastic_dissipation_kpa_m'].append(abs(right['state']['plastic_dissipation_kpa_m']-left['state']['plastic_dissipation_kpa_m']))
    result['fixed_3d_rotation_covariance']=dict(maximum_errors={k:max(v) for k,v in errors.items()},rotation=rotation.tolist(),moving_normal_transport_tested=False)

    (a.output/'rotating.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    summary={k:v for k,v in result.items() if k!='runs'};summary['runs']={k:{n:x for n,x in v.items() if n!='history'} for k,v in result['runs'].items()}
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('rotating runs',len(result['runs']),flush=True)

if __name__=='__main__':main()
