"""Bounded original-open MCC scope control; never soften a numerical stop."""
from pathlib import Path
import argparse
from c03diag_study import read,write,invoke,sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cli',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir()
    case=read(a.source/'stages.json');next(s for s in case['stages'] if s['name']=='open')['target_factor']=1.
    write(out/'stages.json',case);write(out/'ground.json',read(a.source/'ground.json'))
    result=invoke([a.cli.resolve(),'solve-stages',a.matrix.resolve()/'mesh.json',out/'ground.json',out/'stages.json',out/'cli-run.json'],out/'cli.log',900)
    trace=invoke([a.matrix.resolve()/'probe',a.matrix.resolve()/'mesh.json',out/'ground.json',out/'stages.json',out/'raw','--diagnose'],out/'trace.log',900)
    status=read(out/'cli-run.json') if (out/'cli-run.json').exists() else {}
    write(out/'scope-extra.json',dict(cli=result,trace=trace,complete=status.get('complete'),failure_reason=status.get('failure_reason'),
          cli_sha256=sha(a.cli),passed=False,physics_or_threshold_changed=False))

if __name__=='__main__':main()
