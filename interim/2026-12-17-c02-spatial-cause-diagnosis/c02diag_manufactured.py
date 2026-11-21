import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from c02diag_transfer import transfer,apply,norm
import math
import numpy as np


def rule(n):
    x,w=np.polynomial.legendre.leggauss(n)
    x=(x+1)/2;w=w/2
    return np.array([[u,(1-u)*v] for u in x for v in x]),np.array([wu*wv*(1-u) for u,wu in zip(x,w) for v,wv in zip(x,w)])


def children(v):
    a,b,c=v;ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
    return [np.array(t) for t in [(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)]]


def manufactured():
    rs,w=rule(4);targets,tw=rule(12)
    cells=[np.array([[0.,0.],[1.,0.],[0.,1.]])]
    result={m:[] for m in ('nearest','p2_l2','cubic')}
    operators={m:transfer(rs,w,targets,m) for m in result}
    for h in range(4):
        errors={m:[] for m in result}
        weights=[]
        for v in cells:
            x=v[0]+rs@(v[1:]-v[0]);y=v[0]+targets@(v[1:]-v[0])
            exact=y[:,0]**4+y[:,1]**4
            data=x[:,0]**4+x[:,1]**4
            for m,operator in operators.items():
                errors[m].extend((apply(operator,data)-exact).tolist())
            weights.extend((tw*4.**(-h)).tolist())
        for m in result:result[m].append(norm(errors[m],weights))
        if h<3:cells=[c for v in cells for c in children(v)]
    return {m:{'errors':v,'orders':[math.log2(a/b) for a,b in zip(v,v[1:])]} for m,v in result.items()}

if __name__=='__main__':
    import json
    print(json.dumps(manufactured(),indent=2))
