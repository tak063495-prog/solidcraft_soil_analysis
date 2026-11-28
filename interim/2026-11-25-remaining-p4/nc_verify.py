"""Reassemble both ends using the saved stresses and actual P3 mesh."""
import gzip
import json
from pathlib import Path
import sys
import zipfile
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'2026-11-25-remaining-p3'))
from nc_audit import Section, physical_check, prepare

inputs=Path(sys.argv[1])
m=Section(inputs,None)
mesh=json.loads((inputs/'mesh.json').read_bytes())
with zipfile.ZipFile(HERE.parent/'2026-11-23-nc01-coordinate-polish/evidence.zip') as archive:
    run=json.loads(archive.read('run.json'))
cp=run['last_committed']
q,off=prepare(mesh,run['integration_points'])
fixed=np.setdiff1d(np.arange(len(m.base)),m.free)
result={}
for folder in map(Path,sys.argv[2:]):
    summary=json.loads((folder/'summary.json').read_bytes())
    checks={}
    for name,path in [('start',HERE.parents[1]/summary['source']),('last',folder/'last.json.gz')]:
        state=json.load(gzip.open(path,'rt'))
        material=list(cp['material_states'])
        for index,s in zip(m.point_ids,state['active_material_states'],strict=True):
            material[index]=s
        np.testing.assert_array_equal(np.ravel(state['displacement_m'])[fixed],np.ravel(cp['displacement_m'])[fixed])
        force,y,p,pc=physical_check(mesh,run,cp,q,off,material)
        load=m.load(state['load'])
        ratio=np.linalg.norm((load-force.ravel())[m.free])/max(1.,np.linalg.norm(load))
        assert ratio<=m.settings['relative_force_tolerance']
        checks[name]=dict(load=state['load'],force_ratio=ratio,min_p_kpa=p,min_pc_kpa=pc,max_normalized_yield=y)
    result[folder.name]=checks
print(json.dumps(result,indent=2,allow_nan=False))
