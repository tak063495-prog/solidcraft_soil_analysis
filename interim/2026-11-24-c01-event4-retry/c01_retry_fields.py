"""Lossless full-field streams at exact accepted binary load coordinates."""

from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
METHODS = HERE.parent / "2026-11-23-c01-method-improvements"
sys.path.insert(0, str(METHODS))
from run_case import STAGES, context

CONTACT = ("open", "stick", "slip")
RETURNS = (
    "elastic",
    "shear_face",
    "compression_edge",
    "extension_edge",
    "shear_apex",
    "tension_face",
    "tension_edge",
    "tension_apex",
    "mixed",
)
CATEGORIES = {"contact_status", "material_return"}


class Fields:
    def __init__(self, folder, n, side):
        study, mesh, ground, case, _ = context(n)
        self.n = n
        self.side = side
        self.shapes = dict(
            study.FIELDS,
            equivalent_plastic_strain=(324,),
            plastic_work=(324,),
            contact_status=(32,),
            material_return=(324,),
        )
        self.folder = folder
        self.study = study
        self.op = study.prior.operators(mesh, ground, case, "cut") if side == "native" else None
        self.count = 0
        self.digests = {k: hashlib.sha256() for k in self.shapes}
        self.stack = ExitStack()

    def __enter__(self):
        self.folder.mkdir(parents=True)
        try:
            self.streams = {
                k: self.stack.enter_context((self.folder / (k + ".bin")).open("xb"))
                for k in self.shapes
            }
        except BaseException:
            self.stack.close()
            raise
        return self

    def __exit__(self, kind, exc, tb):
        self.stack.close()
        if kind is None:
            files = {
                k: dict(
                    shape=shape,
                    dtype="|u1" if k in CATEGORIES else "<f8",
                    sha256=self.digests[k].hexdigest(),
                )
                for k, shape in self.shapes.items()
            }
            (self.folder / "manifest.json").write_text(
                json.dumps(dict(n=self.n, count=self.count, side=self.side, fields=files), indent=2)
                + "\n",
                encoding="utf-8",
            )

    def write(self, stage, time, values):
        phase, index = divmod(self.count, self.n)
        if (stage, time) != (phase, (index + 1) / self.n) or phase >= len(STAGES):
            raise ValueError("missing or non-accepted prescribed coordinate")
        for name, shape in self.shapes.items():
            arr = np.asarray(values[name])
            if arr.shape != shape or not np.isfinite(arr).all():
                raise ValueError(f"invalid {name} field")
            raw = np.ascontiguousarray(arr, dtype="u1" if name in CATEGORIES else "<f8").tobytes()
            self.streams[name].write(raw)
            self.digests[name].update(raw)
        self.count += 1

    def native(self, stage, time, state):
        # Same extraction and hydrostatic reconstruction as the original C-01 study.
        values = self.study.prior.native_fields(state, self.op)
        material = state["material_states"]
        values.update(
            contact_status=[CONTACT.index(q["status"]) for q in state["joints"]["points"]],
            material_return=[RETURNS.index(q["last_return"]) for q in material],
            equivalent_plastic_strain=[q["equivalent_plastic_strain"] for q in material],
            plastic_work=[q["plastic_work_kpa"] for q in material],
        )
        self.write(stage, time, values)

    def independent(self, name, time, values, status):
        # This reference checks the elastic MC branch in field_state at each sampled GP.
        zeros = np.zeros(self.shapes["material_return"])
        values.update(
            contact_status=status,
            material_return=zeros,
            equivalent_plastic_strain=zeros,
            plastic_work=zeros,
        )
        self.write(STAGES.index(name), time, values)
