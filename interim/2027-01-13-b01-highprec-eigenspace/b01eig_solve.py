"""Offline eigendecomposition of the byte-identified 23c7b0e assemblies."""
import argparse
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from zipfile import ZipFile

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.with_name("2027-01-01-b01-curved-tri6-cloud-followup")
INDEX_BLOB = "f737c7805ec84cccffe817815cfc208099f7a1b4"
INPUT_SHA = "7321ec76510f1b14c1d3654108fe5b602b48fb43d7aca75bb01423405c950583"
INPUT_BLOB = "f8576031d35fd96a3dabc5c3afc35485ac51a49c"
DPS = 45


def blob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def load_inputs():
    index = (OLD / "archive-index.json").read_bytes()
    if blob(index) != INDEX_BLOB:
        raise ValueError("23c7b0e index blob mismatch")
    index = json.loads(index)
    path = OLD / "evidence.zip.part000"
    data = path.read_bytes()
    if (blob(data) != INPUT_BLOB or hashlib.sha256(data).hexdigest() != INPUT_SHA
            or INPUT_SHA != index["archive_sha256"] or len(data) != index["archive_bytes"]):
        raise ValueError("23c7b0e archive identity mismatch")
    cases = []
    from io import BytesIO
    with ZipFile(BytesIO(data)) as archive:
        if sorted(archive.namelist()) != sorted(index["entries"]):
            raise ValueError("incorrect archive member list")
        for name, expected in sorted(index["entries"].items()):
            entry = archive.read(name)
            if len(entry) != expected["size"] or hashlib.sha256(entry).hexdigest() != expected["sha256"]:
                raise ValueError(name)
            if name.startswith("b01follow-formal/assemblies/"):
                case = json.loads(entry)
                if np.shape(case["stiffness"]) != (150, 150):
                    raise ValueError(name)
                cases.append(case)
    if len(cases) != 32 or len({c["id"] for c in cases}) != 32:
        raise ValueError("assembly count")
    return cases


def solve(job):
    import mpmath as mp
    case, output = job
    mp.mp.dps = DPS
    k = np.array(case["stiffness"], dtype=np.float64)
    lam, v = np.linalg.eigh((k + k.T) * 0.5)
    exact = mp.matrix(case["stiffness"])
    sym = (exact + exact.T) / 2
    eigenvalues, vectors = mp.eigsy(sym)
    # mpmath eigsy returns an ascending spectrum with eigenvectors in columns.
    # Store more digits than the working precision for later residual evaluation.
    def fmt(value):
        return mp.nstr(value, DPS + 4)
    record = dict(id=case["id"], shape=case["shape"], rotation=case["rotation"],
                  binary64_values=lam.tolist(), binary64_vectors=v.tolist(),
                  mp_dps=DPS, mp_values=[fmt(x) for x in eigenvalues],
                  mp_vectors=[[fmt(vectors[i,j]) for j in range(150)] for i in range(150)])
    path = output / (case["id"] + ".json")
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(record, stream, separators=(",", ":"), allow_nan=False)
        stream.write("\n")
    print(case["id"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        os.environ[key] = "1"
    cases = load_inputs()
    args.output.mkdir(parents=True, exist_ok=False)
    with ProcessPoolExecutor(args.workers) as pool:
        list(pool.map(solve, [(case, args.output) for case in cases]))
