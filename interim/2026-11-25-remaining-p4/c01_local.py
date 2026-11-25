"""Use the archived local uniqueness probe on the asymmetric cycle."""
from pathlib import Path
here=Path(__file__).resolve().parent
path=here.parent/'2026-11-25-remaining5-diagnosis/c01_branches.py'
code=path.read_text()
old='case["settings"]["attempt_budget_policy"] = "mesh_proportional"'
assert code.count(old)==1
code=code.replace(old,old+'''
for stage in case["stages"]:
    if stage["name"] in ("reverse","reload"):
        stage["target_factor"]={"reverse":-0.075,"reload":0.060}[stage["name"]]
''')
# N1024 already has the requested slip budget. Check one frozen trial at
# the first closure and one at the middle of the actual reload interval.
lo=code.index('    elif name == "reverse":')
hi=code.index('    else:\n        return',lo)
code=code[:lo]+'''    elif name == "recontact":
        if not np.any((v["status"]==0)&(v["new_status"]!=0)):
            return
    elif name in ("reverse","reload"):
        if not v["previous"] <= .5 <= v["progress"]:
            return
'''+code[hi:]
exec(compile(code,str(path),'exec'))
