import json,pathlib
root=pathlib.Path.cwd();c=root/'work/b03cause-final-controls';rows=[]
for path in sorted((c/'cli').glob('*/mesh.json')):
 mesh=json.loads(path.read_text());stage=json.loads((path.parent/'stages.json').read_text());cfg=stage['beam_system']
 rows.append(dict(case=path.parent.name,soil_nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),tri6=len(cfg['plates']),structural_nodes=len(cfg['nodes']),uncondensed_structural_dofs=6*len(cfg['nodes'])))
result=dict(scope='new scientific CLI controls; software-suite historical cases are separately inventoried',internal_free_patches=dict(maximum_tri6=2,maximum_structural_nodes=12,maximum_structural_dofs=72),production_cli_controls=dict(maximum_tet10=max(r['tet10'] for r in rows),maximum_tri6=max(r['tri6'] for r in rows),maximum_soil_nodes=max(r['soil_nodes'] for r in rows),maximum_structural_nodes=max(r['structural_nodes'] for r in rows)),cli_inputs=rows,summary_maximum_keys_scope='maximum_tri6 and maximum_structural_dofs in measurements.json refer only to the internal free patches',new_solves=0)
out=root/'interim/2026-12-26-b03-offset-kink-cause-diagnosis/validation/resource-audit.json';out.write_text(json.dumps(result,indent=2)+'\n');print(result['production_cli_controls'])
