"""Package this reference's result as a proposed diagnostic packet, not a body exporter."""
import argparse
import csv
import json
from pathlib import Path
from compare import packet,validate
from reference import write

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('result',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    data=validate(packet(json.loads(a.result.read_text(encoding='utf-8'))))
    if a.output.exists():raise FileExistsError(a.output)
    if a.output.suffix=='.json':write(a.output,data)
    elif a.output.suffix=='.csv':
        sidecar=a.output.with_suffix('.metadata.json')
        if sidecar.exists():raise FileExistsError(sidecar)
        write(sidecar,data['metadata'])
        with a.output.open('x',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data['rows'][0]));writer.writeheader()
            for row in data['rows']:writer.writerow(dict(row,hydraulic_open=str(row['hydraulic_open']).lower()))
    else:raise ValueError('output .json/.csv required')
