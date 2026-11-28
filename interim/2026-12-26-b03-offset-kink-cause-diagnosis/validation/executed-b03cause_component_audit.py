import json,pathlib,sys
sys.path.insert(0,'interim/2026-12-26-b03-offset-kink-cause-diagnosis')
from b03cause_math import norm
root=pathlib.Path.cwd();base=root/'work/b03cause-final-patches';rows=[]
for p in sorted((base/'assemblies').glob('*.json')):
 a=json.loads(p.read_text());b=json.loads((base/'independent-assemblies'/p.name).read_text());u=[];r=[]
 for i,(x,y) in enumerate(zip(a['global_displacement'],b['global_displacement'])):
  (u if i%6<3 else r).append(x-y)
 rows.append(dict(case=p.stem,translation_difference_l2_m=norm(u),rotation_difference_l2_rad=norm(r),maximum_translation_difference_m=max(map(abs,u)),maximum_rotation_difference_rad=max(map(abs,r))))
result=dict(diagnostic_only=True,passed=False,b03_complete=False,definition='physical components separated from saved native/input-only independent global solutions; mixed-DOF maximum in earlier summary is an algebraic indicator only',cases=rows,maximum_translation_difference_m=max(r['maximum_translation_difference_m'] for r in rows),maximum_rotation_difference_rad=max(r['maximum_rotation_difference_rad'] for r in rows),new_solves=0,new_acceptance_gate=None)
out=root/'interim/2026-12-26-b03-offset-kink-cause-diagnosis/validation';out.mkdir(exist_ok=True);(out/'component-differences.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if k!='cases'})
