from pathlib import Path
import importlib,importlib.metadata,json,site,sys
output=Path(sys.argv[1]);prefix=Path(sys.prefix).resolve();inventory=Path(sys.argv[2])
expected=json.loads(inventory.read_text('utf-8'))['wheels'];versions={v['name']:v['version'] for v in expected.values()}
modules={'gmsh':'gmsh','markdown-it-py':'markdown_it','mdurl':'mdurl','meshio':'meshio','numpy':'numpy','pygments':'pygments','rich':'rich','scipy':'scipy','ijson':'ijson'}
rows={}
for name,version in versions.items():
 assert importlib.metadata.version(name)==version
 module=importlib.import_module(modules[name.lower()]);path=Path(module.__file__).resolve();assert path.is_relative_to(prefix),(name,path,prefix)
 rows[name]=dict(version=version,path=str(path))
assert sys.prefix!=sys.base_prefix and not site.ENABLE_USER_SITE
assert all(Path(path).resolve().is_relative_to(prefix) for path in site.getsitepackages())
with output.open('x',encoding='utf-8') as f:json.dump(dict(passed=True,python=sys.version,prefix=str(prefix),
 same_host_only=True,all9_package_imports_from_new_env=True,user_site_disabled=True,packages=rows,new_PC_or_GPU_acceptance=False),f,indent=2)
print('PASS all9 exact versions imported only from new independent env')
