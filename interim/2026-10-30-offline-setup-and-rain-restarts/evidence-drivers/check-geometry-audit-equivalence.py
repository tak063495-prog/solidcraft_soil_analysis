"""Prove the codemod and compare every original14909 element geometry bitwise."""
from pathlib import Path
import ast, datetime, hashlib, json, os, sys, time

p = Path(__file__).resolve().parent
out = p / 'geometry-audit-equivalence-v1'
out.mkdir()
base = p / 'source635-offline-setup-v2/runtime'
candidate = p / 'geometry-audit-candidate-v1/runtime'
read = lambda q: json.loads(q.read_text('utf-8-sig'))

def sha(q):
    with q.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def dump(q, value):
    with q.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)

class RestoreStandalone(ast.NodeTransformer):
    def visit_Assign(self, node):
        if (len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == 'geometry' and isinstance(node.value, ast.Call)
                and ast.unparse(node.value) == 'ref.Geometry(mesh)'):
            return None
        return self.generic_visit(node)

    def visit_Call(self, node):
        node = self.generic_visit(node)
        if isinstance(node.func, ast.Attribute) and ast.unparse(node.func) == 'geometry.volume':
            node.func = ast.Attribute(value=ast.Name(id='ref', ctx=ast.Load()), attr='volume', ctx=ast.Load())
            node.args.insert(0, ast.Name(id='mesh', ctx=ast.Load()))
        return node

bindings = {str(Path(__file__)): sha(Path(__file__))}
for name in ['tools/verify_projected_stage_package.py', 'tools/verify_high_order.py']:
    original = ast.parse((base / name).read_text('utf-8'))
    updated = ast.parse((candidate / name).read_text('utf-8'))
    # The existing modules contained no other local Geometry use before this
    # codemod, so the normalization proves every other operator/assertion exact.
    normalized = RestoreStandalone().visit(updated)
    assert ast.dump(original) == ast.dump(normalized), name
    for q in [base / name, candidate / name]:
        bindings[str(q)] = sha(q)
mesh_path = p / 'mesh14909/workflow/analysis/ground/mesh.json'
assert read(p / 'mesh14909/study-terminal.json')['passed']
bindings[str(mesh_path)] = sha(mesh_path)
for name in ['tools/fem_reference.py', 'tools/mixed_mcc_io.py', 'tests/test_reference_geometry.py']:
    assert sha(base / name) == sha(candidate / name)
    bindings[str(candidate / name)] = sha(candidate / name)
dump(out / 'driver-child.json', dict(pid=os.getpid(), executable=sys._base_executable,
    script=str(Path(__file__).resolve()), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
for key in ['OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS']:
    os.environ[key] = '1'
sys.path.insert(0, str(candidate / 'tools'))
import numpy as np
import fem_reference as ref
mesh = read(mesh_path)
geometry = ref.Geometry(mesh)
values = 0
start = time.monotonic()
for element in range(len(mesh['tetrahedra'])):
    if (p / 'pause-requested').exists():
        raise SystemExit(2)
    original = ref.volume(mesh, element)
    updated = geometry.volume(element)
    assert set(original) == set(updated)
    for key, a in original.items():
        b = updated[key]
        assert a.dtype == b.dtype and a.shape == b.shape
        assert a.tobytes() == b.tobytes(), (element, key)
        values += a.size
    if element % 1024 == 0:
        print(json.dumps(dict(checked_elements=element+1)), flush=True)
assert all(sha(Path(q)) == h for q, h in bindings.items())
dump(out / 'terminal.json', dict(passed=True, elements=len(mesh['tetrahedra']),
    integration_points=27*len(mesh['tetrahedra']), compared_array_values=values,
    all_shapes_weights_positions_natural_coordinates_B_and_connectivity_bitwise_equal=True,
    both_entire_changed_modules_AST_identical_after_only_local_Geometry_codemod_normalization=True,
    original_rules_constitutive_equations_thresholds_assertions_and_complete_audit_unchanged=True,
    original_sources_and_completed_mesh_unchanged=True, elapsed_seconds=time.monotonic()-start,
    interleaved_validation_not_an_exclusive_performance_benchmark=True,
    full_actual8061_and14909_complete_audits_and_memory_measurement_still_separate=True,
    p6_complete=False, sha256=bindings))
print(json.dumps(dict(passed=True, elements=len(mesh['tetrahedra']), array_values=values)), flush=True)
