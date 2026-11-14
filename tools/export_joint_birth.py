"""Write completed contact-BIRTH result/15 as VTU/PVD/CSV and generation JSON."""
import argparse
from pathlib import Path
from joint_birth_io import export

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh',type=Path)
    parser.add_argument('result',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    export(args.mesh,args.result,args.output)
