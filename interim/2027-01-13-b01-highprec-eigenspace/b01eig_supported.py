"""Resolve the signed worst membrane response into supported spectral subspaces."""
import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import mpmath as mp
import numpy as np

from b01eig_solve import OLD, DPS, load_inputs
from b01eig_measure import D, PREC, array, getcontext, norm, product


def decompose(case, old_solution):
    ids = case["free_dofs"]
    k = array(case["stiffness"])[np.ix_(ids,ids)]
    s = (k+k.T)/2
    vals,vecs = mp.eigsy(mp.matrix([[str(x) for x in row] for row in s]))
    v = np.array([[D(mp.nstr(vecs[i,j],DPS+4)) for j in range(len(ids))]
                  for i in range(len(ids))])
    normal = next(row for row in case["solves"] if row["name"]=="membrane")
    d = array(normal["displacement_rotation"])
    hp = array(old_solution["displacement_rotation"])
    err = (d-hp)[ids].reshape(-1,1)
    coeff = product(v.T,err)
    # Nine free nodes each have one normal rotational coordinate in this plane.
    # Record both ends of the proposed block and every modal weight.
    lower = product(v[:,:9],coeff[:9])
    upper = product(v[:,9:],coeff[9:])
    full_lower, full_upper = 0*d,0*d
    full_lower[ids], full_upper[ids] = lower[:,0],upper[:,0]
    axis = array(case["rotation_matrix"])[:,2]
    normal_basis = np.zeros((len(ids),9),dtype=object)
    for i in range(9):
        normal_basis[6*i+3:6*i+6,i]=axis
    weights = product(v.T,normal_basis)
    backward = product(k,err)[:,0].reshape(-1,6)
    report = dict(id=case["id"], free_dofs=ids, eigenvalues=[mp.nstr(x,DPS+4) for x in vals],
                  normal_coordinate_total_weights=[str(norm(row)**2) for row in weights],
                  modal_error_squared=[str(x*x) for x in coeff[:,0]],
                  projected_error_low9=[str(x) for x in lower[:,0]],
                  projected_error_others=[str(x) for x in upper[:,0]],
                  low9_fraction_error_squared=str(norm(lower)**2/norm(err)**2),
                  reconstruction_error_l2=str(norm(lower+upper-err)),
                  Kff_error_force_l2_kn=str(norm(backward[:,:3])),
                  Kff_error_moment_l2_kn_m=str(norm(backward[:,3:])))
    q = array(case["rotation_matrix"])
    arrays = [product(t.reshape(-1,3),q).reshape(-1,6) for t in (d,hp,full_lower,full_upper)]
    return report,arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    getcontext().prec = PREC
    mp.mp.dps = DPS
    cases = {c["id"]:c for c in load_inputs()}
    records,series = [],[]
    with ZipFile(OLD/"evidence.zip.part000") as archive:
        for rotation in ("identity","rx37"):
            name = "interior-0.06--"+rotation
            precise = json.loads(archive.read(f"b01follow-precision/{name}--membrane--80.json"))
            row,arrays = decompose(cases[name],precise)
            records.append(row)
            series.append(arrays)
    delta = [b-a for a,b in zip(*series,strict=True)]
    node,axis = np.unravel_index(abs(delta[0][:,3:]).argmax(),(25,3))
    components = [x[node,axis+3] for x in delta]
    out = dict(cases=records, node_zero_based=int(node), rotation_axis_xyz=int(axis),
               ordinary_delta_rad=str(components[0]), saved_decimal80_delta_rad=str(components[1]),
               low9_error_delta_rad=str(components[2]), other45_error_delta_rad=str(components[3]),
               signed_closure_rad=str(components[0]-sum(components[1:])))
    with args.output.open("x",encoding="utf-8",newline="\n") as stream:
        json.dump(out,stream,indent=2)
        stream.write("\n")


if __name__=="__main__":
    main()
