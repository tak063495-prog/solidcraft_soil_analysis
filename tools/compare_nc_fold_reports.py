"""Tabulate actual fold evidence; missing series never become a positive verdict."""
import argparse
import json
from pathlib import Path


def compare(paths):
    summaries=[]
    for path in paths:
        report=json.loads(path.read_text(encoding='utf-8'));points=report['points']
        if not points:
            summaries.append(dict(path=str(path),verdict=report['verdict'],accepted_points=0));continue
        d=[p['diagnostic'] for p in points]
        summaries.append(dict(path=str(path),verdict=report['verdict'],accepted_points=len(points),
            peak_lambda=report['peak_lambda'],minimum_sigma_k=min(x['stiffness']['minimum'] for x in d),
            minimum_sigma_augmented=min(x['augmented']['minimum'] for x in d),
            maximum_augmented_condition=(max(x['augmented']['condition_number'] for x in d) if all(x['augmented']['condition_number'] is not None for x in d) else None),
            minimum_normalized_transversality=min(x['normalized_transversality'] for x in d),
            dlambda_ds_range=[min(x['dlambda_ds'] for x in d),max(x['dlambda_ds'] for x in d)],
            maximum_force_residual=max(p['relative_force_residual'] for p in points),
            maximum_force_block_linear_residual=max(p['force_block_linear_residual'] for p in points),
            maximum_augmented_linear_residual=max(p['augmented_linear_residual'] for p in points),
            source_binding=report['source_binding']))
    pairs=[dict(first=a['path'],second=b['path'],peak_lambda_difference=abs(a['peak_lambda']-b['peak_lambda']))
           for i,a in enumerate(summaries) for b in summaries[i+1:] if 'peak_lambda' in a and 'peak_lambda' in b]
    return dict(conclusion_ja='未確定',passed=False,nc01_complete=False,paths=summaries,pairwise_peak_differences=pairs,
                note='Review fixed prehistories, original proportional 1/128 vs 1/256 and input-only reference refinement; no new peak-location tolerance is invented.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reports',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=compare(a.reports)
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
