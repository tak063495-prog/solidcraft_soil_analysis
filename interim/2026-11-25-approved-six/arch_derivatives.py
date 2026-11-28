"""Check the two new element derivatives at a perturbed, unloaded configuration."""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"2026-11-22-impl-t7-physical-instability"))
from curved import CurvedArch, AnalyticChordArch


def check(model, second):
    n = len(model.free)
    x = np.arange(n)
    z = np.r_[1e-5*np.sin(0.73*x), 0.]
    direction = np.cos(0.421*x)
    direction /= np.linalg.norm(direction)
    step = 2e-7
    dz = np.r_[direction,0]*step
    left, right = model({},z-dz), model({},z+dz)
    if second:
        approximate = ((right.k-left.k)@direction)/(2*step)
        exact = model.tangent_action_derivative(z,direction)
    else:
        approximate = (right.residual-left.residual)/(2*step)
        exact = model({},z).k@direction
    relative = np.linalg.norm(exact-approximate)/np.linalg.norm(exact)
    assert relative < 1e-6
    return float(relative)


if __name__ == "__main__":
    print(json.dumps(dict(curved_12_force_derivative=check(CurvedArch(12),False),
                          chord_96_tangent_derivative=check(AnalyticChordArch(96),True)),indent=2))
