"""Re-encode already saved fields losslessly; do not invoke any solver."""
import argparse
import copy
import io
import json
from pathlib import Path
import shutil
import numpy as np
import study


def compact(source, output):
    output.mkdir(parents=True, exist_ok=False)
    for name in ('mesh.json', 'ground.json', 'stages.json', 'timings.json'):
        shutil.copyfile(source/name, output/name)
    totals = dict(chunks=0, accepted_states=0, original_bytes=0, compact_bytes=0,
                  every_field_bitwise_roundtrip_verified=True, solver_rerun=False)
    for g, local in study.ROUTES:
        key = f'g{g}-l{local}'
        src = source/key/'accepted'; dst = output/key/'accepted'
        dst.mkdir(parents=True)
        shutil.copyfile(source/key/'result.json', output/key/'result.json')
        old = json.loads((src/'manifest.json').read_text())
        shutil.copyfile(src/'manifest.json', output/key/'original-manifest.json')
        manifest = copy.deepcopy(old); manifest['encoding'] = study.ENCODING
        for chunk in manifest['chunks']:
            raw = (src/chunk['file']).read_bytes()
            assert study.sha(raw) == chunk['sha256']
            with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
                fields = {k: study.decode(archive[k], old['encoding']) for k in archive.files}
            assert set(fields) == set(study.FIELDS)
            with (dst/chunk['file']).open('xb') as stream:
                np.savez_compressed(stream, **{k: study.encode(v) for k, v in fields.items()})
            encoded = (dst/chunk['file']).read_bytes()
            with np.load(io.BytesIO(encoded), allow_pickle=False) as archive:
                for name, original in fields.items():
                    np.testing.assert_array_equal(study.decode(archive[name]).view('<u8'), original.view('<u8'))
            totals['chunks'] += 1; totals['accepted_states'] += chunk['count']
            totals['original_bytes'] += len(raw); totals['compact_bytes'] += len(encoded)
            chunk['sha256'] = study.sha(encoded); chunk['bytes'] = len(encoded)
        study.write(dst/'manifest.json', manifest)
        print('losslessly compacted', key, totals, flush=True)
    study.write(output/'compaction.json', totals)
    result = study.measure(output)
    study.write(output/'measurements.json', result)
    return totals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    compact(args.source.resolve(), args.output.resolve())
