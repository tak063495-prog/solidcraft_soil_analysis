"""Input-only closed forms. SI m/s, kN, kPa; z upward, compression positive.

No production modules, subprocesses, solver states or result files are used.
Derivation and verified repository equation sources: README.ja.md.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Column:
    top: float = 3.0
    interface: float = 1.2  # pre-existing/new body boundary
    water: float = 2.1
    gamma_dry: float = 17.0
    gamma_sat: float = 20.0
    gamma_w: float = 9.81
    surcharge: float = 5.0
    k0: float = 0.55  # prescribed, NOT inferred from friction angle/OCR

    def pressure(self, z):
        return self.gamma_w * np.maximum(self.water - np.asarray(z), 0.0)

    def weight(self, z):
        return np.where(np.asarray(z) < self.water, self.gamma_sat, self.gamma_dry)

    def total_vertical(self, z):
        z = np.asarray(z)
        wet = np.maximum(np.minimum(self.top, self.water) - z, 0.0)
        dry = self.top - z - wet
        top_water = self.gamma_w * max(self.water - self.top, 0.0)
        return self.surcharge + top_water + self.gamma_dry * dry + self.gamma_sat * wet

    def fields(self, z, pressure=None):
        z = np.asarray(z)
        p = self.pressure(z) if pressure is None else np.asarray(pressure)
        sv = self.total_vertical(z)
        ev = sv - p
        eh = self.k0 * ev
        return dict(pressure=p, total_vertical=sv, effective_vertical=ev,
                    effective_horizontal=eh, total_horizontal=eh+p)


@dataclass(frozen=True)
class Darcy:
    edges: tuple = (0.0, 1.0, 3.0)
    conductivity: tuple = (1e-5, 4e-5)
    bottom_head: float = 4.4
    top_head: float = 3.6
    gamma_w: float = 9.81

    def __post_init__(self):
        if (len(self.edges) != len(self.conductivity)+1 or
                not np.all(np.diff(self.edges) > 0) or
                not np.all(np.asarray(self.conductivity) > 0)):
            raise ValueError('strictly increasing layers and positive conductivity required')

    @property
    def flux(self):
        return (self.bottom_head-self.top_head)/sum(
            (b-a)/k for a, b, k in zip(self.edges[:-1], self.edges[1:], self.conductivity))

    def head(self, z):
        z = np.asarray(z)
        resistance = np.zeros_like(z, dtype=float)
        for a, b, k in zip(self.edges[:-1], self.edges[1:], self.conductivity):
            resistance += np.clip(z-a, 0, b-a)/k
        return self.bottom_head-self.flux*resistance

    def k(self, z):
        i = np.clip(np.searchsorted(self.edges, z, side='right')-1, 0, len(self.conductivity)-1)
        return np.asarray(self.conductivity)[i]

    def pressure(self, z):
        p = self.gamma_w*(self.head(z)-np.asarray(z))
        if np.any(p < 0):
            raise ValueError('negative saturated pressure: no clipping/free-surface inference')
        return p

    def pressure_gradient(self, z):
        return self.gamma_w*(-self.flux/self.k(z)-1.0)

    def seepage(self, z):
        return self.gamma_w*self.flux/self.k(z)  # upward force density kN/m3


def cases():
    return {
        'fully_submerged': (Column(water=3.8), None),
        'partial_new_layer': (Column(water=2.1), None),
        'water_at_birth_interface': (Column(water=1.2), None),
        'dry_new_layer': (Column(water=0.6), None),
        'darcy_upward': (Column(water=3.6), Darcy()),
        'darcy_downward': (Column(water=4.4), Darcy(bottom_head=3.6, top_head=4.4)),
        'darcy_three_layers': (Column(water=3.6), Darcy(
            edges=(0., .7, 1.8, 3.), conductivity=(1e-5, 6e-5, 2e-5))),
    }
