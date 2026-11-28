"""Sequential continuous-face inclusions using a polynomial center and past-history bounds."""
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P11 = HERE.parent / "2026-11-26-remaining-p11"
sys.path.insert(0, str(P11))
from face_matrix import Reduction, base

EPS = np.finfo(float).eps


def norm_hi(x):
    return np.linalg.norm(x, axis=-1)*(1+32*EPS)+1e-95


def setup(sign):
    folder = base.prev.P8/"c02"/sign
    blocks, law = base.prev.geometry(folder)
    refined = []
    for ids, _, _, areas, tiles, sgn in blocks:
        children=np.array([c for tri in tiles for c in base.prev.p6.split(tri)])
        w=base.prev.transform(children)
        trace2=2*np.sum(np.max(abs(w),axis=1)**2,axis=1)*(1+256*EPS)
        refined.append((ids,w.reshape(-1,15),np.max(abs(w).sum(axis=2)),sgn,np.repeat(areas/4,4),trace2))
    area=np.concatenate([r[-2] for r in refined])
    assert len(area)==6080 and abs(area.sum()-1)<1e-12
    reducers=[]
    for weight in ("1","4"):
        os.environ["P11_TEST_WEIGHT"]=weight
        reducers.append(Reduction(refined,area))
    data=np.load(folder/"kinematics.npz")
    return refined, reducers, np.sqrt(np.concatenate([r[-1] for r in refined]))*(1+32*EPS), law, data


def enclosure(refined,u):
    jumps=[]
    for ids,w,amp,sgn,_,_ in refined:
        left,right=u.reshape(-1,3)[ids]
        jump=(right-left)*[1,1,sgn]
        coeff=(w@jump).reshape(-1,15,3)
        pad=192*EPS*amp*(abs(left).max()+abs(right).max())+1e-95
        jumps.append((coeff.min(1)-pad,coeff.max(1)+pad))
    return tuple(map(np.concatenate,zip(*jumps)))


def expand(bounds, width):
    l,h=bounds
    return (np.nextafter(l-width[:,None],-np.inf),np.nextafter(h+width[:,None],np.inf))


def defects(previous,current,D,width,first,kn,kt,mu):
    """Bound |tau-kt*j_t| for every possible previous history and material partition."""
    a,b=previous
    c,f=current
    closed=(b[:,2]<=0)&(f[:,2]<=0)
    q=kt*norm_hi(np.maximum(-c[:,:2],f[:,:2]))
    disklo=mu*kn*np.maximum(-f[:,2],0)*(1-32*EPS)
    diskhi=mu*kn*np.maximum(-c[:,2],0)*(1+32*EPS)
    elastic=closed & ((q+D)*(1+64*EPS)<disklo)
    open_=c[:,2]>0
    # Zero history and proportional loading: P at the endpoint for every partition.
    # On the closed side f_cl=1; the open side has exactly zero final traction.
    if first:
        return np.maximum(0,q-disklo)*(1+64*EPS),elastic,open_
    delta=np.maximum(b-c,f-a)*(1+32*EPS)+1e-90
    count=np.maximum(1,np.ceil(norm_hi(delta)/width*(1+32*EPS)))
    # Compare each projected update to kt*j_t. Projection is nonexpansive;
    # discarded translation has magnitude <= |delta_t| (1-f_cl) over the path.
    # The disk excess at every substep is bounded by max |kt*j_t| minus min radius.
    maxq=kt*norm_hi(np.maximum.reduce([-a[:,:2],b[:,:2],-c[:,:2],f[:,:2]]))
    min_disk=mu*kn*np.maximum(-np.maximum(b[:,2],f[:,2]),0)*(1-32*EPS)
    fraction_lower=base.interval_fraction(a[:,2],b[:,2],c[:,2],f[:,2])[0]
    path=D+kt*norm_hi(delta[:,:2])*(1-fraction_lower)+count*np.maximum(0,maxq-min_disk)
    bound=np.minimum(q+diskhi,path)*(1+128*EPS)+1e-90
    # For a wholly closed segment the old trial is admissible at its start.
    # If the extrapolated end is admissible, cone convexity covers all subdivisions.
    return np.where(elastic,D,np.where(open_,q,bound)),elastic,open_


def integrate(red,u,t,phase,closed_center,gb_hi,d_nom,hist,kn,kt):
    g=red.g
    bulk=g.k@u+g.initial-g.load-((1-t)*(g.initial-g.load) if phase==0 else 0.)
    force=bulk.reshape(-1,3).copy()
    center_err=np.zeros_like(force)
    past=np.zeros_like(force)
    scale=float(np.linalg.norm(abs(g.k)@abs(u)+2*(abs(g.initial)+abs(g.load))))
    for sl,_,ids,_,mass,_,amax in red.blocks:
        # Exact degree-eight weak-form integral of the polynomial spring field.
        vals=(closed_center[sl].astype(float)@mass).reshape(15,15)
        j=u.reshape(-1,3)[ids[1]]-u.reshape(-1,3)[ids[0]]
        spring=vals@(j*[kt,kt,kn])
        scale+=float(np.linalg.norm(abs(vals)@(abs(j)*[kt,kt,kn])))
        force[ids[0]]-=spring
        force[ids[1]]+=spring
        dt=d_nom[sl]*closed_center[sl]
        err=amax.T@np.column_stack((dt,dt,kn*np.maximum(gb_hi[sl],0)*closed_center[sl]))
        h=amax.T@hist[sl]
        for side in ids:
            center_err[side]+=err
            past[side,:2]+=h[:,None]
    free=red.free
    guard=4096*EPS*(1+scale)
    center=float(np.linalg.norm((abs(force)+center_err).ravel()[free])*(1+1e-10)+guard)
    past_upper=float(np.linalg.norm(past.ravel()[free])*(1+1e-10))
    return center,past_upper,guard


def first_lower(red,eligible,opened,kn,kt,mu):
    # Closed-side projected disk derivative has a PSD tangential part and
    # cross derivative bounded by mu*kn. Completing the square against w*kn
    # spends at most mu^2*kn/(4*w) on both tangential coordinates. The zero
    # open branch satisfies this same bound, including at one-sided switches.
    losses=np.where(eligible,kt,np.where(opened,0.,-mu**2*kn/(4*red.weight)*(1+1e-12)))
    normal=np.where(eligible,red.weight*kn,0.)
    M=np.zeros((75,75))
    for sl,pos,_,_,mass,_,_ in red.blocks:
        for j,values in enumerate((losses,losses,normal)):
            M[np.ix_(3*pos+j,3*pos+j)]+=(values[sl]@mass).reshape(15,15)
    margins={}
    for floor,G in red.floors.items():
        margin=float(np.linalg.eigvalsh(G+(M+M.T)/2)[0]-(.1+1e-8*max(kn,kt)*red.weight))
        margins[floor]=margin
        if margin>0:
            return floor/red.weight,margins
    return 0.,margins


def run(sign,out):
    refined,reds,jump,law,data=setup(sign)
    kn,kt,mu=[law[k] for k in ("normal_stiffness_kpa_m","tangential_stiffness_kpa_m","friction_coefficient")]
    clock,U=data["clock"],data["u"]
    zero=np.zeros((len(jump),3))
    prev_center=prev_ball=(zero,zero)
    dnom=dtrue=np.zeros(len(jump))
    prev_radius=0.
    rows=[]
    area=np.concatenate([r[-2] for r in refined])
    closed_translation=np.zeros_like(jump)
    anchor_bounds=[]
    folder=base.prev.P8/"c02"/sign
    slip=np.load(folder/"slip.npz")
    radii=json.loads((HERE/"design.json").read_text())["C02"]["nested_test_radii_m"]
    old=[[json.loads(s) for s in (prefix/f"{sign}-1e-08.jsonl").read_text().splitlines()]
         for prefix in (P11/"c02",P11/"c02/weighted")]
    for i,((phase,t,width),u) in enumerate(zip(clock,U,strict=True)):
        first=i==0
        center=enclosure(refined,u)
        dc,_,opened=defects(prev_center,center,dnom,width,first,kn,kt,mu)
        actual,_,_=defects(prev_ball,center,dtrue,width,first,kn,kt,mu)
        # At a fixed current displacement normal traction has no history term.
        # Both tangential values have the SAME kt*j_t reference here.
        hist=np.zeros_like(dc) if first else np.where(opened,0.,(dc+actual)*(1+32*EPS))
        fc,fh,guard=integrate(reds[0],u,t,phase,~opened,center[1][:,2],dc,hist,kn,kt)
        residual=(fc+fh)*(1+16*EPS)
        tests=[]
        carried=None
        for radius in radii:
            ball=expand(center,jump*radius)
            _,elastic,op=defects(prev_ball,ball,dtrue,width,first,kn,kt,mu)
            margins={}
            if first:
                floors=[]
                for red in reds:
                    v,m=first_lower(red,elastic,op,kn,kt,mu)
                    margins[str(red.weight)]=m
                    floors.append(v)
                lower=max(floors)
            elif (elastic|op).all():
                # On every child the tangent is constant PSD on this ball.
                lower=next(iter(reds[0].floors))  # K - 700 I certified in Reduction
            else:
                assert all(path[i]["phase"]==phase and path[i]["t"]==t for path in old)
                lower=max(path[i]["lower"] or 0 for path in old)
            ok=residual<lower*radius
            tests.append(dict(radius=radius,lower_kn_m=lower,elastic=int(elastic.sum()),open=int(op.sum()),
                              one_sided_cells=(dict(closed=int((ball[1][:,2]<=0).sum()),open=int(op.sum()),
                                  both=int(((ball[0][:,2]<=0)&(ball[1][:,2]>0)).sum())) if first else None),
                              schur_margins=margins,included=bool(ok)))
            if ok:
                carried=(residual/lower)*(1+64*EPS)
                break
        rows.append(dict(step=i+1,phase=int(phase),t=t,previous_radius_m=prev_radius,
                         initial_exact_zero=first,center_upper_kn=fc,past_upper_kn=fh,assembly_guard_kn=guard,
                         total_upper_kn=residual,tests=tests,certified_radius_m=carried))
        if carried is None:
            break
        prev_ball_new=expand(center,jump*carried)
        dtrue,_,_=defects(prev_ball,prev_ball_new,dtrue,width,first,kn,kt,mu)
        delta=np.maximum(prev_ball[1][:,:2]-prev_ball_new[0][:,:2],prev_ball_new[1][:,:2]-prev_ball[0][:,:2])
        closed_translation=(closed_translation+norm_hi(delta))*(1+8*EPS)
        if t*64==int(t*64):
            # |slip| <= sum of closed translation lengths + |tau|/kt.
            norm_bound=closed_translation+mu*kn/kt*np.maximum(-prev_ball_new[0][:,2],0)
            active=norm_hi(slip["slip"][len(anchor_bounds)])
            anchor_bounds.append(float((area@norm_bound+abs(slip["weight"]).ravel()@active)*(1+1e-10)))
        dnom=dc
        prev_center,prev_ball,prev_radius=center,prev_ball_new,carried
    chain=sum(r["certified_radius_m"] is not None for r in rows)
    with (out/(sign+".jsonl")).open("x") as stream:
        for row in rows:
            print(json.dumps(row,allow_nan=False),file=stream)
    ratio=None
    if chain==len(clock):
        assert len(anchor_bounds)==448
        ref=json.loads((folder.parent/"comparison.json").read_text())[sign]["bands"]["all"]
        ratio=float(np.linalg.norm(anchor_bounds)*(1+1e-10)/ref["own_active_signed_integral"]["denominator"])
    (out/(sign+"-summary.json")).write_text(json.dumps(dict(consecutive=chain,total=len(clock),
             threshold=.001,reequilibrated_full_history_relative_upper=ratio,first=rows[0],stop=rows[-1]),indent=2)+"\n")


if __name__=="__main__":
    run(sys.argv[1],Path(sys.argv[2]))
