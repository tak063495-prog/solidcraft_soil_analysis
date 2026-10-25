"""Cached-start midpoint integration for the independent MCC reference.

Only a derivative evaluated at the identical committed state and active set is
reused during event bisection. The constitutive law, event brackets, accepted
updates and drift limits are unchanged. No solution state is cached across steps.
"""
import numpy as np


def advance_cached(self,state,displacement,force_rate,count,progress=None):
    terms=self.constitutive_terms
    t=0.;state=state.copy();u=displacement.copy();plastic=np.zeros(len(state),dtype=bool)
    maximum_drift=0.;events=0;accepted=0;minimum_step=1.;step=1/count
    while t<1-1e-13:
        plastic=self.select(state,force_rate,plastic)
        _,start_rate,_,_=self.rate(state,force_rate,plastic)
        def midpoint(h):
            du,dy,_,_=self.rate(state+start_rate*h/2,force_rate,plastic)
            return state+dy*h,du*h
        h=min(step,1-t);candidate,du=midpoint(h)
        yf=terms(candidate,self.parameters)[0]
        _,_,potential,_=self.rate(candidate,force_rate,plastic)
        crossed=(~plastic)&(yf>1e-9);unloaded=plastic&(potential<-1e-11);next_plastic=plastic.copy()
        if np.any(crossed)|np.any(unloaded):
            lo=0.;hi=h
            for _ in range(28):
                mid=(lo+hi)/2;trial,_=midpoint(mid)
                fy=terms(trial,self.parameters)[0]
                _,_,ld,_=self.rate(trial,force_rate,plastic)
                if np.any((~plastic)&(fy>1e-9)) or np.any(plastic&(ld<-1e-11)):hi=mid
                else:lo=mid
            h=hi;candidate,du=midpoint(h)
            fy=terms(candidate,self.parameters)[0];_,_,ld,_=self.rate(candidate,force_rate,plastic)
            next_plastic|=(~plastic)&(fy>=0.999e-9);next_plastic&=ld>=-0.999e-11
            events+=1
        actual_yield=terms(candidate,self.parameters)[0]
        drift=float(np.max(np.maximum(actual_yield,0.)))
        if np.any(next_plastic):drift=max(drift,float(np.max(np.abs(actual_yield[next_plastic]))))
        if drift>2e-6:
            step=h/2
            if step<1e-9:raise ValueError('independent yield drift cannot be controlled without a projection')
            continue
        if h<1e-11:raise ValueError('independent event integration stalled')
        state=candidate;u+=du;plastic=next_plastic;t+=h;accepted+=1;minimum_step=min(minimum_step,h);maximum_drift=max(maximum_drift,drift)
        if progress is not None and accepted%100==0:progress(t,accepted,events)
    return state,u,dict(accepted_steps=accepted,events=events,maximum_normalized_yield_drift=maximum_drift,
                       minimum_progress_step=minimum_step,rate_assemblies=self.calls,
                       independent_pcg_iterations=self.linear_iterations,dense_fallbacks=self.dense_fallbacks)
