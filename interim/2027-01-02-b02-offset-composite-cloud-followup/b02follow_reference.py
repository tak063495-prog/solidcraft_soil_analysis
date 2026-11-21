"""Only the independent plate quadrature changes; no production/soil/beam override."""
import argparse
from pathlib import Path
from b02follow_common import quadrature, legacy

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--order',type=int,choices=(4,5,6),required=True)
    a=p.parse_args()
    with quadrature(a.order):
        legacy['reference'].solve(a.folder,a.output)
