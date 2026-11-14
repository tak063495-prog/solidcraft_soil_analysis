"""Re-audit saved small-cycle fields; never run a solver or certify histories."""
import argparse
import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
from verify_joint_water import audit
from verify_wet_structures import read, sha, without_time
from verify_wet_joint_cycle import path_evidence
from wet_joint_cycle_reference import compare, differences, native_fields


def remeasure(evidence):
    evidence = Path(evidence)
    manifest = read(Path(__file__).with_name('archive-sha256.json'))
    for name, expected in manifest['entries_sha256'].items():
        assert sha(evidence/name) == expected, ('evidence SHA mismatch', name)
    old = read(evidence/'measurements.json'); measured = copy.deepcopy(old)
    assert old['passed'] is False and old['reference_uses_native_results'] is False
    counts = old['reference_global_and_local_subincrements']
    refs = [read(evidence/f'reference-{n}.json') for n in counts]
    order7 = read(evidence/'reference-order7.json')
    natives = []
    for name in ['original']+[f'h{n}' for n in counts]:
        folder = evidence/name; mesh = read(folder/'mesh.json')
        ground = fem.normalize_ground(read(folder/'ground.json'))
        native = read(folder/'stage-results.json'); natives.append(native)
        assert native['provenance']['solver_binary_sha256'] == old['solver_sha256']
        measured['native_routes'][name] = path_evidence(native)
        try:
            measured['native_vs_reference'][name] = compare(native, refs[-1])
            measured['existing_displacement_and_traction_gates'][name] = True
        except AssertionError:
            measured['existing_displacement_and_traction_gates'][name] = False
            expected = copy.deepcopy(refs[-1]); keys = native_fields(native)['frames'][0].keys()
            for f in expected['frames']:
                for key in list(f.keys()-keys):
                    del f[key]
            measured['native_vs_reference'][name] = differences(native_fields(native), expected)
        measured['native_stored_field_and_output_audits'][name] = audit(mesh, ground, native, folder/'results')
        exported = read(folder/'results/export-verification.json')
        assert exported['result_sha256'] == sha(folder/'stage-results.json')
        assert exported['mesh_sha256'] == sha(folder/'mesh.json')
    measured['independent_step_sensitivity'] = differences(*refs)
    measured['independent_volume_order_sensitivity'] = differences(refs[-1], order7)
    measured['native_step_sensitivity'] = differences(native_fields(natives[1]), native_fields(natives[2]))
    for row in old['restart_checks']:
        name = row['stage']; folder = evidence/'original'
        assert sha(folder/f'checkpoint-{name}.json') == row['checkpoint_sha256']
        assert without_time(read(folder/f'resumed-{name}.json')) == without_time(natives[0])
    assert measured == old, 'saved measurements do not reproduce exactly'
    return measured


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True); args = p.parse_args()
    result = remeasure(args.evidence)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    print('Saved measurements reproduced exactly; passed=false retained.')
