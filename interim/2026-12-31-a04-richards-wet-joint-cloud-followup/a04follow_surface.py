"""Complete Q through the independent dedicated hydrostatic ray quadrature."""
import argparse
from contextlib import contextmanager
import shutil
from a04follow_common import *
from a04follow_joint_order import binding, controls, layout

@contextmanager
def ray_order():
    import hydrostatic_reference as hydro
    original=hydro.rule
    def rule(n):
        x,w=original(n)
        return x[::-1].copy(),w[::-1].copy()
    hydro.rule=rule
    try:yield
    finally:hydro.rule=original

def worker(folder):
    s,r,_=binding('hydro')
    with controls(True),ray_order():
        model=r.Model(*[read(folder/(k+'.json')) for k in ('mesh','ground','case')])
        write(folder/'a04follow-independent-layout.json',layout(model))
        s.independent(folder,4)
    import hydrostatic_reference as hydro
    write(folder/'ray-point-set-audit.json',dict(rules={str(n):dict(count=n,point_weight_pairs=sorted(zip(*[v.tolist() for v in hydro.rule(n)])),reversal_keeps_exact_pairs=True) for n in (4,12,24)},
        fully_submerged_planar_faces=True,reversed_inner_and_outer=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('--formal',type=Path);p.add_argument('--output',type=Path);p.add_argument('--worker',type=Path);a=p.parse_args()
    if a.worker:worker(a.worker);return
    a.output.mkdir(parents=True,exist_ok=False)
    sources={p.name:sha(p) for p in HERE.glob('*.py')};write(a.output/'source-binding.json',sources)
    snap=a.output/'source-snapshot';snap.mkdir()
    for name in sources:shutil.copy2(HERE/name,snap/name)
    paths=[]
    for e,f in ((0,0),(0,1),(1,0),(1,1)):
        name=f'e{e}f{f}q1';folder=a.output/name;folder.mkdir()
        for key in ('mesh','ground','case'):shutil.copy2(a.formal/'order/hydro'/name/(key+'.json'),folder/(key+'.json'))
        proc=invoke([sys.executable,Path(__file__),'--worker',folder],folder,'independent')
        paths.append(dict(path=name,process=proc))
    assert sources=={p.name:sha(p) for p in HERE.glob('*.py')}
    write(a.output/'result.json',dict(passed=False,a04_complete=False,diagnostic_only=True,execution_complete=True,fresh_paths=4,paths=paths))

if __name__=='__main__':main()
