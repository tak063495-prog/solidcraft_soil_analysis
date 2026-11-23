"""Compare complete accepted load grids or their measured common prefix, field by field."""

import hashlib
import json
from math import prod

import numpy as np

from c01_retry_fields import CATEGORIES, STAGES


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


class History:
    def __init__(self, folder, side, n, stages):
        self.folder = folder
        self.meta = json.loads((folder / "manifest.json").read_bytes())
        self.n = n
        self.count = self.meta["count"]
        expected = sum(int((s["last_accepted"] or 0) * n) for s in stages)
        if (self.meta["side"], self.meta["n"], self.count) != (side, n, expected):
            raise ValueError("field count differs from accepted stage records")
        for name, row in self.meta["fields"].items():
            path = folder / (name + ".bin")
            dtype = "|u1" if name in CATEGORIES else "<f8"
            if (
                row["dtype"] != dtype
                or path.stat().st_size != self.count * prod(row["shape"]) * np.dtype(dtype).itemsize
                or sha256(path) != row["sha256"]
            ):
                raise ValueError(f"invalid field stream {name}")

    def grid(self, name, n):
        if self.n % n:
            raise ValueError("non-nested load grids")
        row = self.meta["fields"][name]
        if not self.count:
            return np.empty((0, *row["shape"]), dtype=row["dtype"])
        values = np.memmap(
            self.folder / (name + ".bin"),
            dtype=row["dtype"],
            mode="r",
            shape=(self.count, *row["shape"]),
        )
        stride = self.n // n
        return values[stride - 1 :: stride]


def compare_pair(a, b, n, category=False):
    if a.shape != b.shape:
        raise ValueError("field shape or accepted coordinates differ")
    if category:
        return dict(mismatches=int(np.count_nonzero(a != b)), components=a.size)
    if not len(a):
        return dict(
            absolute_l2=None,
            reference_l2=None,
            relative_l2=None,
            maximum_absolute=None,
            maximum_location=None,
            components=0,
        )
    d = a - b
    norm = float(np.linalg.norm(d))
    ref = float(np.linalg.norm(b))
    ix = tuple(map(int, np.unravel_index(np.argmax(np.abs(d)), d.shape)))
    stage, index = divmod(ix[0], n)
    return dict(
        absolute_l2=norm,
        reference_l2=ref,
        relative_l2=norm / max(ref, 1e-12),
        maximum_absolute=abs(float(d[ix])),
        maximum_location=[STAGES[stage], (index + 1) / n, *ix[1:]],
        components=d.size,
    )


def observed(a, b, c):
    if not len(a):
        return dict(
            coarse_norm=None, fine_norm=None, observed_order=None, cosine=None, trend="unobserved"
        )
    ab, bc = a - b, b - c
    x, y = float(np.linalg.norm(ab)), float(np.linalg.norm(bc))
    p = float(np.log2(x / y)) if min(x, y) > 0 else None
    cosine = float(np.vdot(ab, bc) / (x * y)) if min(x, y) > 0 else None
    return dict(
        coarse_norm=x,
        fine_norm=y,
        observed_order=p,
        cosine=cosine,
        trend=(
            "identical"
            if x == y == 0
            else "undefined"
            if p is None
            else "decreasing_aligned"
            if p > 0 and cosine > 0
            else "withheld"
        ),
    )


def align(series, field, n):
    grids = [s.grid(field, n) for s in series]
    count = min(map(len, grids))
    return [g[:count] for g in grids]


def compare(histories):
    fields = sorted(next(iter(histories.values())).meta["fields"])
    if any(sorted(s.meta["fields"]) != fields for s in histories.values()):
        raise ValueError("field inventories differ")
    refinements, solvers = {}, {}
    for mode in ("local", "ws", "local-ws", "local-active"):
        for side in ("native", "independent"):
            series = [histories[side, n, mode] for n in (256, 512, 1024)]
            common = min(s.count // (s.n // 256) for s in series)
            entry = dict(common_times=common, complete=common == len(STAGES) * 256, fields={})
            for field in fields:
                cat = field in CATEGORIES
                pairs = {}
                for i, n in enumerate((256, 512)):
                    a, b = align(series[i : i + 2], field, n)
                    pairs[str(n)] = dict(times=len(a), **compare_pair(a, b, n, cat))
                a, b, c = align(series, field, 256)
                metrics = (
                    dict(
                        coarse_mismatches=int(np.count_nonzero(a != b)),
                        fine_mismatches=int(np.count_nonzero(b != c)),
                        components=a.size,
                    )
                    if cat
                    else observed(a, b, c)
                )
                entry["fields"][field] = dict(adjacent=pairs, common=metrics)
                if field == "plastic_slip":
                    phases = {}
                    for i, stage in enumerate(STAGES):
                        window = slice(i * 256, (i + 1) * 256)
                        if len(a[window]):
                            phases[stage] = dict(
                                times=len(a[window]), **observed(a[window], b[window], c[window])
                            )
                    entry["fields"][field]["by_stage"] = phases
            refinements[f"{side}-{mode}"] = entry
        for n in (256, 512, 1024):
            series = [histories[side, n, mode] for side in ("native", "independent")]
            common = min(s.count for s in series)
            solvers[f"{n}-{mode}"] = dict(
                times=common,
                complete=common == len(STAGES) * n,
                fields={
                    field: compare_pair(*align(series, field, n), n, field in CATEGORIES)
                    for field in fields
                },
            )
    return dict(refinements=refinements, solvers=solvers)
