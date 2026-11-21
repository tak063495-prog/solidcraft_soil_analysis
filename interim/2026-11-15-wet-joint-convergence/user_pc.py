"""Explicit user-PC extension. Never invoked by CI or the cloud study."""
import argparse
from pathlib import Path
import study
from study import read,write,refine,fem,ROOT
from measure import compare,spatial
from joint_water_reference import Water


def run(output,solver,levels,steps,maximum,native_timeout):
    output.mkdir(parents=True,exist_ok=False)
    mesh=read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    ground=fem.normalize_ground(read(ROOT/'examples/ground-joint-water-mc.json'))
    case=fem.normalize_stage_case(read(ROOT/'examples/stages-joint-water-cycle.json'))
    write(output/'ground.json',ground);write(output/'stages.json',case)
    for h in range(max(levels)+1):
        write(output/f'mesh-h{h}.json',mesh)
        if h<max(levels):mesh,audit=refine(mesh);write(output/f'refinement-h{h+1}.json',audit)
    routes=[dict(id=f'h{h}-g{g}-m{maximum:g}-cut',h=h,g=g,local=1,rule='cut',maximum=maximum,native=True,user_pc=True,native_timeout_seconds=native_timeout) for h in levels for g in steps]
    write(output/'routes.json',routes)
    solver=study.freeze_solver(output,solver)
    failures={}
    for r in routes:
        try:study.solve_route(output,r,solver)
        except Exception as e:
            failures[r['id']]=dict(passed=False,complete=False,error=repr(e));write(output/r['id']/'failure.json',failures[r['id']])
    comparisons={}
    def path(h,g,k):return output/f'h{h}-g{g}-m{maximum:g}-cut'/k
    for kind in ('independent','native'):
        for h in levels:
            for a,b in zip(steps,steps[1:]):
                if (path(h,a,kind)/'manifest.json').exists() and (path(h,b,kind)/'manifest.json').exists():
                    comparisons[f'{kind}/h{h}/g{a}__{b}']=compare(path(h,a,kind),path(h,b,kind))
        for a,b in zip(levels,levels[1:]):
            assert b==a+1
            systems=[Water(read(output/f'mesh-h{h}.json'),ground,case).contact for h in (a,b)]
            for g in steps:
                if (path(a,g,kind)/'manifest.json').exists() and (path(b,g,kind)/'manifest.json').exists():
                    comparisons[f'{kind}/h{a}__{b}/g{g}']=spatial(path(a,g,kind),path(b,g,kind),*systems)
    write(output/'measurements.json',dict(passed=False,diagnostic_completed=not failures,failures=failures,comparisons=comparisons,
                                        full_history_accepted=False,reason='User-PC extension measurements only; unchanged existing criteria and unresolved pointwise spatial mapping.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--levels',type=int,nargs='+',default=[0,1,2]);p.add_argument('--steps',type=int,nargs='+',default=[64,128,256]);p.add_argument('--maximum',type=float,default=5e-5);p.add_argument('--native-timeout-seconds',type=int,default=14400)
    a=p.parse_args();assert a.levels==sorted(set(a.levels)) and a.steps==sorted(set(a.steps)) and min(a.levels)>=0
    assert a.maximum>0 and __import__('math').isfinite(a.maximum) and a.native_timeout_seconds>0
    run(a.output.resolve(),a.solver.resolve(),a.levels,a.steps,a.maximum,a.native_timeout_seconds)
