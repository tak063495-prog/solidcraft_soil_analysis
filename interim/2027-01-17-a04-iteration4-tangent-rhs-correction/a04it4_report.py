"""Compare actual/replayed systems, mapped exactly to ascending physical free DOFs."""
import argparse
from decimal import localcontext
import math
from pathlib import Path
import sys

from a04it4_decimal import D, decimal_array, elastic, np, solve
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'2026-12-31-a04-richards-wet-joint-cloud-followup'))
from a04follow_common import read,rows,write


def peak(x):
    return float(max(map(abs,x.flat)))


def metric(x):
    idx = np.unravel_index(np.argmax(abs(x)),x.shape)
    return dict(max_abs=float(abs(x[idx])),index=list(map(int,idx)),signed=str(x[idx]))


def canonical(basis, record):
    # Signed coordinate selection only. No floating projection or SVD in comparisons.
    assert np.all(np.isin(basis,[-1,0,1])) and np.all((basis!=0).sum(axis=0)==1)
    ids = np.flatnonzero(np.any(basis!=0,axis=1))
    assert len(ids)==basis.shape[1]
    perm = np.argmax(abs(basis[ids]),axis=1)
    assert len(set(perm))==len(ids)
    sign = basis[ids,perm].astype(int)
    a = decimal_array(record['matrix'])[np.ix_(perm,perm)]*sign[:,None]*sign
    b = decimal_array(record['rhs'])[perm]*sign
    x = decimal_array(record['dx'])[perm]*sign
    return ids,a,b,x


def side(a,b,x,K,R):
    high = solve(a,b)
    defect = b-a@x
    scale = max(abs(a).sum(axis=1))*max(abs(x))+max(abs(b))
    return dict(tangent_vs_elastic=metric(a-K),rhs_vs_elastic=metric(b-R),
        linear_error_m=metric(x-high),linear_defect_kn=metric(defect),
        backward_error_inf=float(max(abs(defect))/scale),
        tangent_error_in_peak_matrix_ulps=peak(a-K)/math.ulp(peak(K)),
        rhs_error_in_peak_rhs_ulps=peak(b-R)/math.ulp(peak(R)),
        linear_error_in_peak_solution_ulps=peak(x-high)/math.ulp(peak(high))),high


def report(root,out,digits):
    out.mkdir(parents=True,exist_ok=False)
    folders = [root/'000',root/'010']
    with localcontext() as ctx:
        ctx.prec = digits
        K,seed = elastic(read(folders[0]/'mesh.json'),read(folders[0]/'ground.json'))
        write(out/'elastic.json',dict(matrix=K.astype(str),seed=seed.astype(str)))
        result,paths = {},[]
        for folder in folders:
            history = list(rows(folder/'native.jsonl'))
            it = next(r for r in history if r['kind']=='iterate' and r['i']==4)
            accepted = next(r for r in history if r['kind']=='accepted')
            linear = next(r for r in history if r['kind']=='linear')
            indep = read(folder/'independent-fixed.json')
            t = np.zeros((len(linear['rows']),len(linear['rhs'])))
            for g,terms in enumerate(linear['rows']):
                for j,c in terms:
                    t[g,j] = c
            ids,a,b,x = canonical(t,linear)
            ids2,c,d,y = canonical(np.array(indep['basis']),indep)
            assert np.array_equal(ids,ids2) and set(indep['native_status'])=={'open'}
            assert not np.any(indep['contact_tangent']) and len(indep['gap'])==32
            u = decimal_array(it['u']).ravel()
            load = decimal_array(it['target'])
            k = K[np.ix_(ids,ids)]
            R = (load-K@u-seed)[ids]
            n,xn = side(a,b,x,k,R)
            i,yi = side(c,d,y,k,R)
            # Signed telescoping terms in the same coordinates for x_native-x_independent.
            cross = solve(a,d)
            tangent,rhs,linear_err = cross-yi,xn-cross,(x-xn)-(y-yi)
            assert accepted['alpha']==1
            next_u = decimal_array(accepted['u']).ravel()
            ix_next = decimal_array(indep['full_step'])
            fifth = next(r for r in history if r['kind']=='iterate' and r['i']==5)
            assert np.array_equal(fifth['u'],accepted['u'])
            last = [r for r in history if r['kind']=='iterate'][-1]
            def reaction(row):
                return decimal_array(row['internal'])-decimal_array(row['target'])
            paths.append((a,b,x,xn,u,next_u,ids,load,reaction(fifth),reaction(last)))
            own = read(folder/'independent-path.json')
            end = [r for r in history if r['kind']=='iterate'][-1]
            result[folder.name] = dict(free_dofs=len(ids),native=n,independent_fixed=i,
                peak_exact_tangent=peak(k),peak_exact_rhs=peak(R),peak_correction=peak(xn),
                native_minus_fixed=dict(tangent=metric(a-c),rhs=metric(b-d),dx=metric(x-y),
                    tangent_term_m=metric(tangent),rhs_term_m=metric(rhs),solver_term_m=metric(linear_err),
                    closure_m=peak(x-y-tangent-rhs-linear_err)),
                native_update_rounding_m=metric(next_u[ids]-u[ids]-x),
                independent_trial_rounding_m=metric(ix_next[ids]-u[ids]-y),
                native_alpha=accepted['alpha'],native_iteration5_relative=accepted['relative'],
                native_stop=end['i'],independent_stop=own[-1]['i'],
                independent_relative=[r['relative'] for r in own],
                states_open=len(indep['gap']),minimum_gap_m=min(indep['gap']))
            write(out/(folder.name+'.json'),{key:value.astype(str) for key,value in
                  dict(ids=ids,K=k,rhs=R,native_exact_solve=xn,independent_exact_solve=yi,
                       cross_solve=cross,tangent_term=tangent,rhs_term=rhs,linear_term=linear_err).items()})
        a,b,x,xh,u,v,ids,L,f0,end0 = paths[0]
        c,d,y,yh,w,z,ids2,M,f1,end1 = paths[1]
        assert np.array_equal(ids,ids2)
        swap = solve(c,b)
        # 010 minus 000, change matrix first with b000 fixed, then change RHS.
        before = (w-u)[ids]
        matrix_step,rhs_step,linear_step = swap-xh,yh-swap,(y-yh)-(x-xh)
        update = (z-v)[ids]-before-(y-x)
        delta = z-v
        net = (K@delta)-(M-L)
        ideal0 = u.copy()
        ideal1 = w.copy()
        ideal0[ids] += solve(K[np.ix_(ids,ids)], (L-K@u-seed)[ids])
        ideal1[ids] += solve(K[np.ix_(ids,ids)], (M-K@w-seed)[ids])
        result['010_minus_000'] = dict(tangent=metric(c-a),rhs=metric(d-b),
            incoming_u_m=metric(before),incoming_load_kn=metric(M-L),dx_m=metric(y-x),
            tangent_term_m=metric(matrix_step),rhs_term_m=metric(rhs_step),
            linear_term_m=metric(linear_step),update_term_m=metric(update),
            closure_m=peak(delta[ids]-before-matrix_step-rhs_step-linear_step-update),
            incoming_plus_rhs_term_m=metric(before+rhs_step),
            candidate_u_m=metric(delta),candidate_reaction_kn=metric(net),
            observed_fifth_reaction_kn=metric(f1-f0),fifth_prediction_error_kn=metric(f1-f0-net),
            observed_endpoint_reaction_kn=metric(end1-end0),
            ideal_newton_candidate_m=metric(ideal1-ideal0))
        write(out/'orders.json',{key:value.astype(str) for key,value in
            dict(incoming=before,tangent_term=matrix_step,rhs_term=rhs_step,linear_term=linear_step,
                 update_rounding=update,delta_u=delta,delta_force=net,ideal_newton_delta=ideal1-ideal0).items()})
        write(out/'summary.json',result)


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--digits',type=int,choices=(40,55),default=40)
    a=p.parse_args()
    report(a.raw,a.output,a.digits)
