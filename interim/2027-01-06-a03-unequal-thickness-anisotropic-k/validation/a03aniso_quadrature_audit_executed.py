"""Input-only assembly audit; zero new equilibrium solves or accepted states."""
import sys,math,json,hashlib
from pathlib import Path
root=Path.cwd();here=root/'interim/2027-01-06-a03-unequal-thickness-anisotropic-k';sys.path.insert(0,str(here))
from a03aniso_binding import reference,study,controls
rows=[]
for spec in study.read(here/'plan.json')['specs']:
 if spec['n']!=8:continue
 folder=root/'work/a03aniso/formal'/spec['name']
 with controls(spec['volume'],False,1):model=reference.Model(*[study.read(p) for p in study.paths(folder)])
 counts=[len(op[0]['weight']) for op in model.ops];assert counts==[spec['volume']**3]*12
 rows.append(dict(path=spec['name'],volume_order=spec['volume'],actual_assembly_points_per_element=counts,
                  output_material_points_per_element=[len(q['weight']) for q in model.sample],
                  element_volumes_m3=[math.fsum(map(float,op[0]['weight'])) for op in model.ops]))
(here/'validation/quadrature-contract.json').write_text(json.dumps(dict(assembly_audit_only=True,new_equilibrium_solves=0,new_accepted_states=0,
 legacy_reference_header_is_constant_GL125=True,effective_rule_is_volume_order_and_observed_counts=True,rows=rows),indent=2)+'\n')
print('verified',len(rows),'input assemblies; GL5=125/GL7=343; zero new solves')
