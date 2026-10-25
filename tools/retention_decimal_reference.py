"""Independent 100-digit Decimal retention reference; not a PDE oracle."""
import argparse
from decimal import Decimal as D, localcontext
from pathlib import Path
import hashlib
import json


def reference():
  with localcontext() as c:
    c.prec=100
    tr,ts,alpha,n,ss=map(D,['0.08','0.45','1.2','1.6','0.00001']);m=1-1/n
    rows=['pressure_head_m,saturation,theta,deficit,capacity,relative_conductivity']
    for text in ['-1e20','-1e8','-100','-10','-1','-0.1','-1e-4','-1e-12','-1e-30','0','1','100']:
        p=D(text)
        if p>=0:s=D(1);theta=ts;deficit=D(0);capacity=ss;k=D(1)
        else:
            x=-alpha*p;s=(1+x**n)**(-m);theta=tr+(ts-tr)*s;deficit=(ts-tr)*(1-s)
            capacity=(ts-tr)*m*n*alpha*x**(n-1)*(1+x**n)**(-m-1)
            k=s.sqrt()*(1-(1-s**(1/m))**m)**2
        rows.append(','.join(map(str,[p,s,theta,deficit,capacity,k])))
    return '\n'.join(rows)+'\n'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',type=Path);parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args();data=reference()
    assert args.verify or args.output_dir,'select --verify and/or --output-dir'
    if args.verify:assert args.verify.read_text(encoding='utf-8-sig')==data,'Decimal fixture differs'
    if args.output_dir:
        args.output_dir.mkdir(parents=True,exist_ok=False);target=args.output_dir/'decimal-reference.csv'
        with target.open('x',encoding='utf-8') as stream:stream.write(data)
        files=[Path(__file__).resolve(),target,Path(__file__).resolve().parents[1]/'src/retention.rs']
        with (args.output_dir/'provenance.json').open('x',encoding='utf-8') as stream:json.dump({'precision_decimal_digits':100,'scope':'Retention kernel only; no PDE solution is inferred.',
            'sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}},stream,indent=2)
    print('Independent Decimal retention reference verified')
