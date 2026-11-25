"""Account for corrector calls, including excluded diagnostics and P3 secant replays."""
from pathlib import Path
import json

folder=Path(__file__).resolve().parent/'nc'
counts={str(p.parent.relative_to(folder)):sum(1 for _ in p.open()) for p in folder.rglob('trials.jsonl')}
counts['fine P3 replay including fold']=400
counts['auxiliary P3 replay including fold']=256
counts['partial calls at two interrupts (conservative)']=2
counts['zero-increment tangent assemblies for excluded starts']=3
counts['fixed-point material-only replays (17 at two units)']=34
counts['same-parent audits (two methods at two terminal routes, allow repartition)']=8
total=sum(counts.values())
assert total<=(2065*96+11)//12
print(json.dumps(dict(counts=counts,total_upper=total,mesh_budget=16520),indent=2))
