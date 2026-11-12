"""Prepare larger existing-input jobs; never increase geometry or solver permission."""
import argparse
from pathlib import Path
from b03_cases import make,write

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path)
    p.add_argument('--kind',choices=('kink-zero','kink-split','step-split','normal-sign-pair','height-step'),default='kink-zero')
    p.add_argument('--level',type=int,choices=range(4),default=2);p.add_argument('--steps',type=int,choices=(4,8,16),default=4)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);m,g,c,refinement=make(a.kind,a.level,a.steps)
    for k,v in [('mesh',m),('ground',g),('stages',c),('refinement',refinement)]:write(a.output/(k+'.json'),v)
    write(a.output/'scope.json',dict(passed=False,b03_complete=False,kind=a.kind,nodes=len(m['nodes']),tet10=len(m['tetrahedra']),
        plates=len(c['beam_system']['plates']),independent_direct_reference_within_soil_limit=3*len(m['nodes'])<=6500,
        note='Preparation only. Existing capacity/geometry refusal is retained; no clipping, extra constraints, or independent limit override.'))
