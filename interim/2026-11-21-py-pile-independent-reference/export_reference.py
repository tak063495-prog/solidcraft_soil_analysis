"""Export an identified 1D case, or an explicitly labelled verification case.

This exporter cannot represent the historical square 3D pile. It never infers
spring stiffness from soil E and never calls the application or old reference.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_pile import ContractError, read_json, require, validate_metadata
from pile_reference import axial_exact, beam_exact, evaluate, finite_elements


def verification_metadata(kind):
    axial = kind == "axial"
    fields = ({"u": ("m", "axial_displacement"), "N": ("kN", "compression_force"),
               "shaft_reaction": ("kN/m", "axial_soil_on_pile")} if axial else
              {"w": ("m", "lateral_displacement"), "slope": ("rad", "dw_ds"),
               "M": ("kN*m", "minus_EI_d2w_ds2"), "V": ("kN", "minus_EI_d3w_ds3"),
               "soil_reaction": ("kN/m", "lateral_soil_on_pile")})
    return dict(schema="pile-field-comparison/1", units="m-kN-kPa",
                case=dict(id="verification-only-"+kind, model="1d-axial-linear-springs" if axial else "1d-winkler-eb-free-tip",
                          geometry=dict(length_m=6., area_m2=1., inertia_m4=1/12),
                          materials=dict(pile_young_kpa=30_000_000.,
                                         distributed_k_kn_m2=30_000_000/36 if axial else 4*2_500_000/6**4,
                                         tip_k_kn_m=5_000_000. if axial else 0.),
                          loads=dict(head_load_kn=100. if axial else 10., head_moment_kn_m=0., distributed_load="none"),
                          boundary="axial elastic tip" if axial else "head moment zero; tip M=V=0",
                          interface="independent linear springs; verification parameters, NOT identified 3D mapping",
                          water="dry", reference_state="load increment; no gravity increment; small strain linear"),
                basis=dict(origin_m=[0., 0., 0.], axes=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
                signs=({"u": "positive global -Z (down)", "N": "compression positive; N=-EA*du/ds",
                        "shaft_reaction": "positive global -Z; force of soil on pile=-ka*u"} if axial else
                       {"w": "positive global +X", "slope": "dw/ds; s increases global -Z",
                        "M": "M=-EI*d2w/ds2", "V": "V=-EI*d3w/ds3; V(0)=-H",
                        "soil_reaction": "positive global +X; force of soil on pile=-kh*w"}),
                frame="verification-only full head load relative to unloaded spring model",
                quantities={q: dict(unit=unit, component=definition, location="section at station; FE right limit except tip left",
                                    definition=definition+"; s=8-z increasing down; u positive down; w positive global X")
                            for q, (unit, definition) in fields.items()},
                expected_ids={q: [f"s{i:03d}" for i in range(17)] for q in fields},
                provenance=dict(source="explicit verification-only parameters; README.ja.md",
                                sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


def export(meta, kind, method, n=16):
    validate_metadata(meta)
    expected_model = "1d-axial-linear-springs" if kind == "axial" else "1d-winkler-eb-free-tip"
    require(meta["case"]["model"] == expected_model, "3D/other model cannot use a 1D reference")
    c = meta["case"]
    require(set(c["geometry"]) == {"length_m", "area_m2", "inertia_m4"}, "unsupported geometry")
    require(set(c["materials"]) == {"pile_young_kpa", "distributed_k_kn_m2", "tip_k_kn_m"}, "unsupported constitutive conditions")
    require(set(c["loads"]) == {"head_load_kn", "head_moment_kn_m", "distributed_load"} and
            c["loads"]["head_moment_kn_m"] == 0 and c["loads"]["distributed_load"] == "none", "unsupported loads")
    require(c["interface"] in {"independent linear springs",
                               "independent linear springs; verification parameters, NOT identified 3D mapping"},
            "unsupported interface; no continuum-to-spring conversion")
    require(c["water"] == "dry" and c["reference_state"] == "load increment; no gravity increment; small strain linear",
            "only dry, linear load increments are supported")
    require(c["boundary"] == ("axial elastic tip" if kind == "axial" else "head moment zero; tip M=V=0"),
            "unsupported boundary condition")
    require(meta["basis"] == dict(origin_m=[0., 0., 0.], axes=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
            "exporter supports explicit global XYZ world basis only")
    sample = verification_metadata(kind)
    require(meta["quantities"] == sample["quantities"] and meta["signs"] == sample["signs"] and
            meta["expected_ids"] == sample["expected_ids"], "unsupported sampling/quantity/sign contract")
    g, m, f = c["geometry"], c["materials"], c["loads"]
    length, spring, load, tip = g["length_m"], m["distributed_k_kn_m2"], f["head_load_kn"], m["tip_k_kn_m"]
    require(all(isinstance(x, (int, float)) and not isinstance(x, bool) and x > 0
                for x in [length, g["area_m2"], g["inertia_m4"], m["pile_young_kpa"]]), "positive geometry/material required")
    stiffness = m["pile_young_kpa"] * g["area_m2" if kind == "axial" else "inertia_m4"]
    require(kind == "axial" or tip == 0., "beam tip spring not supported")
    s = np.linspace(0., length, 17)
    if method == "exact":
        values = (axial_exact(s, length, stiffness, spring, load, tip) if kind == "axial" else
                  beam_exact(s, length, stiffness, spring, load))
    else:
        require(method == "fe", "method must be exact or fe")
        solution, _ = finite_elements(kind, n, length, stiffness, spring, load, tip)
        values = evaluate(kind, solution, s, length, stiffness, spring)
    meta = json.loads(json.dumps(meta))
    meta["provenance"] = dict(source=f"standalone {kind} {method}, n={n if method == 'fe' else 'not applicable'}",
                              sha256=hashlib.sha256(Path(__file__).with_name("pile_reference.py").read_bytes()).hexdigest())
    rows = [dict(id=f"s{i:03d}", quantity=q, x_m=5., y_m=5., z_m=float(8-s[i]), value=float(values[i,j]))
            for j,q in enumerate(sample["quantities"]) for i in range(len(s))]
    return dict(metadata=meta, rows=rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--verification-only", action="store_true")
    mode.add_argument("--metadata", type=Path)
    parser.add_argument("--kind", choices=["axial", "beam"], required=True)
    parser.add_argument("--method", choices=["exact", "fe"], required=True)
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        meta = verification_metadata(args.kind) if args.verification_only else read_json(args.metadata)
        result = export(meta, args.kind, args.method, args.n)
        with args.output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    except (ValueError, KeyError, TypeError, OSError) as exc:
        parser.exit(2, f"reference refused: {exc}\n")


if __name__ == "__main__":
    main()
