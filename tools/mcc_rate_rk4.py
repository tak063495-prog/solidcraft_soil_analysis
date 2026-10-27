"""Fourth-order alternative for the independent MCC reference, without projection.

Uses the same analytical rate law, active-set selection and event/drift limits
as the midpoint reference. Every event bracket recomputes its RK4 trial from
the unchanged step start. No constitutive state is corrected onto the ellipse.
"""
import numpy as np


def advance_rk4(self,state,displacement,force_rate,count,progress=None):
    if type(count) is not int or count<2:raise ValueError('at least two RK4 steps required')
    terms=self.constitutive_terms
    t=0.;state=state.copy();u=displacement.copy();plastic=np.zeros(len(state),dtype=bool)
    maximum_drift=0.;events=0;accepted=0;minimum_step=1.;step=1/count
    while t<1-1e-13:
        plastic=self.select(state,force_rate,plastic)
        u1,k1,_,_=self.rate(state,force_rate,plastic)
        def trial(h):
            u2,k2,_,_=self.rate(state+k1*h/2,force_rate,plastic)
            u3,k3,_,_=self.rate(state+k2*h/2,force_rate,plastic)
            u4,k4,_,_=self.rate(state+k3*h,force_rate,plastic)
            return state+(k1+2*k2+2*k3+k4)*(h/6),(u1+2*u2+2*u3+u4)*(h/6)
        h=min(step,1-t);candidate,du=trial(h)
        yf=terms(candidate,self.parameters)[0]
        _,_,potential,_=self.rate(candidate,force_rate,plastic)
        crossed=(~plastic)&(yf>1e-9)&(potential>1e-12);unloaded=plastic&(potential<-1e-11);next_plastic=plastic.copy()
        if np.any(crossed)|np.any(unloaded):
            from scipy.optimize import brentq
            def margin(length):
                local,_=trial(length)
                fy=terms(local,self.parameters)[0]
                _,_,ld,_=self.rate(local,force_rate,plastic)
                # min implements AND, max implements OR; the zero set is
                # the exact existing event predicate with unchanged units/limits.
                score=np.full(len(local),-np.inf)
                score[~plastic]=np.minimum(fy[~plastic]-1e-9,ld[~plastic]-1e-12)
                score[plastic]=-ld[plastic]-1e-11
                return float(np.max(score))
            left=margin(0.);right=margin(h)
            if not left<=0<right:raise ValueError('independent RK4 event bracket is inconsistent')
            h=brentq(margin,0.,h,xtol=min(1e-12,h/2**28),rtol=4*np.finfo(float).eps,maxiter=80)
            candidate,du=trial(h)
            fy=terms(candidate,self.parameters)[0];_,_,ld,_=self.rate(candidate,force_rate,plastic)
            next_plastic|=(~plastic)&(fy>=0.999e-9);next_plastic&=ld>=-0.999e-11
            events+=1
        actual_yield=terms(candidate,self.parameters)[0]
        drift=float(np.max(np.maximum(actual_yield,0.)))
        if np.any(next_plastic):drift=max(drift,float(np.max(np.abs(actual_yield[next_plastic]))))
        if drift>2e-6:
            step=h/2
            if step<1e-9:raise ValueError('independent RK4 yield drift cannot be controlled without a projection')
            continue
        if h<1e-11:raise ValueError('independent RK4 event integration stalled')
        state=candidate;u+=du;plastic=next_plastic;t+=h;accepted+=1;minimum_step=min(minimum_step,h);maximum_drift=max(maximum_drift,drift)
        if progress is not None and accepted%100==0:progress(t,accepted,events)
    return state,u,dict(integration_method='rk4',event_localizer='brentq_continuous_original_margins',accepted_steps=accepted,events=events,
        maximum_normalized_yield_drift=maximum_drift,minimum_progress_step=minimum_step,
        rate_assemblies=self.calls,independent_pcg_iterations=self.linear_iterations,dense_fallbacks=self.dense_fallbacks)
