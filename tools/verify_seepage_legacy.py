"""Check previous MCC/IDH and prescribed-water results on identical inputs."""
import argparse
import json
from pathlib import Path
from verify_checkpoint import invoke
from verify_mcc_regression import verify
from verify_water_legacy import compare


def run(solver,output):
    compare(solver,output);root=Path(__file__).resolve().parents[1]
    for name,folder in [('staged-water-mc','runs/water-acceptance-20261024-081027-d482faa3/cad-mc'),('staged-water-mcc','runs/water-acceptance-20261024-081027-d482faa3/cad-mcc')]:
        baseline=root/folder;target=output/(name+'-result.json')
        inputs=[baseline/'ground/mesh.json',baseline/'ground/ground.json',baseline/'stages.json']
        invoke([str(solver),'solve-stages',*map(str,inputs),str(target)],output/(name+'.log'))
        verify(baseline/'stage-results.json',target,output/(name+'-comparison.json'))
    with (output/'seepage-legacy-verification.json').open('x') as stream:json.dump({'schema':'soil-fem-seepage-legacy/1','passed':True,'exact_previous_cases':6},stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args();run(a.solver.resolve(),a.output_dir.resolve())
