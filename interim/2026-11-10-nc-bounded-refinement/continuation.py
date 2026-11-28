"""Manual bounded NC continuation; no production solver or acceptance change."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
from collections import deque
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS = ROOT/'interim/2026-11-10-nc-mean-stress-release'
sys.path[:0] = [str(ROOT/'tools'), str(PREVIOUS)]
from linear_evidence_codec import encode_rows


class Streams:
    """Save every full field; all linear summaries and the last 128 actual matrices.

    The complete raw-system stream is hashed, not claimed reconstructible. Its
    per-record SHA plus matrix SHA/rhs/solution permit exact restart comparisons.
    """
    def __init__(self):
        self.items = {}
        self.tail = deque(maxlen=128)

    def append(self, path, row):
        import gzip
        raw = (json.dumps(row, allow_nan=False)+'\n').encode()
        linear = path.name == 'linear-systems.jsonl'
        if path not in self.items:
            self.items[path] = dict(sha=hashlib.sha256(), size=0, rows=0,
                file=gzip.open(str(path)+'.gz', 'xb', compresslevel=1))
        item = self.items[path]
        item['sha'].update(raw); item['size'] += len(raw); item['rows'] += 1
        if linear:
            self.tail.append(row)
            matrix = np.asarray(row['matrix'], dtype='<f8')
            rhs, solution = np.asarray(row['rhs']), np.asarray(row['solution'])
            summary = {k: v for k, v in row.items() if k != 'matrix'}
            summary.update(matrix_sha256=hashlib.sha256(matrix.tobytes()).hexdigest(),
                matrix_shape=list(matrix.shape), raw_record_sha256=hashlib.sha256(raw).hexdigest(),
                independently_recomputed_true_relative_residual=float(np.linalg.norm(matrix@solution-rhs)/max(np.linalg.norm(rhs), 1e-30)))
            item['file'].write((json.dumps(summary, allow_nan=False)+'\n').encode())
        else:
            item['file'].write(raw)

    def close(self, out):
        manifest = {}
        for path, item in self.items.items():
            item['file'].close()
            linear = path.name == 'linear-systems.jsonl'
            manifest[str(path.relative_to(out))] = dict(raw_sha256=item['sha'].hexdigest(), raw_size=item['size'],
                rows=item['rows'], encoding='all-summaries-last-128-matrices' if linear else 'gzip-jsonl',
                complete_raw_rows_reconstructible=not linear)
        if self.tail:
            encode_rows(iter(self.tail), out/'terminal-linear-systems.jsonl.f64.gz')
        (out/'streaming-digests.json').write_text(json.dumps(manifest, indent=2)+'\n')


def refined_solve(original, matrix, rhs):
    """Same matrix and float64 gate; extended residual and equilibrated corrections."""
    answer = original(matrix, rhs)
    norm = max(np.linalg.norm(rhs), 1e-30)
    error = float(np.linalg.norm(matrix@answer-rhs)/norm)
    history = [error]
    if error > 1e-12:
        row_scale = np.max(np.abs(matrix), axis=1)
        if np.any(row_scale == 0):
            return answer, history
        first = matrix/row_scale[:, None]
        column_scale = np.max(np.abs(first), axis=0)
        if np.any(column_scale == 0):
            return answer, history
        scaled = first/column_scale[None, :]
        extended = matrix.astype(np.longdouble)
        for _ in range(8):
            residual = np.asarray(rhs.astype(np.longdouble)-extended@answer.astype(np.longdouble), dtype=float)
            correction = original(scaled, residual/row_scale)/column_scale
            candidate = answer+correction
            candidate_error = float(np.linalg.norm(matrix@candidate-rhs)/norm)
            history.append(candidate_error)
            if candidate_error >= error:
                break
            answer, error = candidate, candidate_error
    return answer, history


def run(args):
    module = importlib.import_module('bounded_study' if args.map == 'native' else 'bounded_rate_study')
    zero = importlib.import_module('study')
    # The file is imported by absolute path under a distinct name by reproduce.py.
    # Avoid this driver's basename shadowing the historical lambda-zero module.
    streams = Streams()
    module.append = streams.append
    zero.append = streams.append
    original_write = module.write
    def write(path, value):
        original_write(path, value)
        if path.name == 'restart.json':
            coordinate = value['payload']['coordinate']
            threshold = .97 if args.order == 'complement-first' else .515
            tail = path.parent/'tail.json'
            if threshold <= coordinate < 1. and not tail.exists():
                # Copy the actual sealed checkpoint; never edit/reseal a payload.
                tail.write_bytes(path.read_bytes())
    module.write = write
    original_model = module.Model
    source = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    class Model(original_model):
        def __init__(self, controls, out):
            super().__init__(controls, out)
            self.source['followup'] = dict(driver_sha256=source, linear_policy=args.linear_policy,
                numpy_version=np.__version__, longdouble_mantissa_bits=np.finfo(np.longdouble).nmant,
                observer_source_ref='f19ce642647207a87b9dda6cc31b006b85b5936c')
    module.Model = Model
    original = np.linalg.solve
    counter = 0
    def solve(matrix, rhs):
        nonlocal counter
        counter += 1
        if args.linear_policy == 'refined':
            answer, history = refined_solve(original, matrix, rhs)
        else:
            answer = original(matrix, rhs)
            history = [float(np.linalg.norm(matrix@answer-rhs)/max(np.linalg.norm(rhs), 1e-30))]
        streams.append(args.output/'linear-refinement.jsonl', dict(call=counter, original_relative_residual=history[0],
            candidate_relative_residuals=history[1:], returned_relative_residual=float(np.linalg.norm(matrix@answer-rhs)/max(np.linalg.norm(rhs), 1e-30)),
            matrix_sha256=hashlib.sha256(np.asarray(matrix, dtype='<f8').tobytes()).hexdigest()))
        return answer
    np.linalg.solve = solve
    try:
        module.run(args)
    finally:
        np.linalg.solve = original
        streams.close(args.output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map', choices=['native', 'independent'], required=True)
    p.add_argument('--order', choices=['mean-first', 'complement-first'], default='complement-first')
    for name in ['evidence', 'binary', 'reference', 'resume']:
        p.add_argument('--'+name, type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--step-size', type=float, choices=[1/128, 1/256], default=1/128)
    p.add_argument('--material-steps', type=int, choices=[4, 8], default=4)
    p.add_argument('--minimum-step', type=float, choices=[2**-22, 2**-28], default=2**-28)
    p.add_argument('--count', type=int, default=4096)
    p.add_argument('--budget-seconds', type=float, default=900)
    p.add_argument('--linear-policy', choices=['plain', 'refined'], default='refined')
    a = p.parse_args()
    if not 0 < a.budget_seconds <= 1800 or not 0 < a.count <= 4096:
        p.error('bounded controls required')
    run(a)
