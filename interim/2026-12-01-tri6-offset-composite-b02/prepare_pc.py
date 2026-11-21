"""Prepare existing-input large PC diagnostics; no large independent acceptance."""
import argparse
from pathlib import Path
from cases import inputs, write

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--level',type=int,default=3);p.add_argument('--steps',type=int,default=16)
    p.add_argument('--offset',type=float,default=.05);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    m,g,c,audits=inputs(a.level,a.steps,a.offset,True)
    for name,value in [('mesh',m),('ground',g),('stages',c),('refinement',audits)]:write(a.output/(name+'.json'),value)
    write(a.output/'scope.json',dict(nodes=len(m['nodes']),tet10=len(m['tetrahedra']),plates=len(c['beam_system']['plates']),
          soil_dofs=3*len(m['nodes']),expected_accepted_points=5*a.steps,
          independent_direct_supported=3*len(m['nodes'])<=6500,passed=False,b02_complete=False,
          note='Native-only large diagnostic when direct reference exceeds existing 6500 soil DOF bound. Storage/geometry guards remain unchanged.'))
