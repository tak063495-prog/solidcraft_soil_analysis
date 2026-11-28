"""Re-evaluate the p18 enclosure while retaining its additive leaf supports."""
import importlib.util
import inspect
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P18 = HERE.with_name("2026-11-28-remaining-p18")
spec = importlib.util.spec_from_file_location("p18_face", P18/"face.py")
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)

paired = previous.p16.paired


def measured(red, refined, displacements, radii, leaf_box, kt, km, A, ix, direction, dest):
    source = inspect.getsource(paired)
    def replace(find, replacement):
        nonlocal source
        assert source.count(find) == 1, find
        source = source.replace(find, replacement)

    if dest.stem.startswith("plus-") and int(dest.stem.split("-")[1]) >= 833:
        replace("    norm_a =", "    restrict = (mats[:,None] @ restrict[None]).reshape(64,15,15)\n    norm_a =")
    replace("    for block, mesh in zip(red.blocks, refined, strict=True):",
            "    rows = []\n    for block_id,(block, mesh) in enumerate(zip(red.blocks, refined, strict=True)):")
    replace("            component += size*", "            leaf_component = size*")
    replace("            tiles += len(cells)", """            component += leaf_component
            leaf_force = np.zeros_like(midpoint_force)
            for nodes,side in zip(ids,(-1,1),strict=True):
                leaf_force[nodes,:2] += side*m
            rows.append((block_id,k,size*(amps @ rad),size*(directional+amp_guard*rad).sum(),
                         leaf_component,A@leaf_force.ravel()[red.free]))
            tiles += len(cells)""")
    replace("    return action,", """    if rows:
        block, leaf, scalar, along_leaf, components, mid = map(np.asarray,zip(*rows))
        np.savez_compressed(dest.with_name(dest.stem+"-supports.npz"),
            block_leaf=np.column_stack((block,leaf)),scalar_m=scalar,directional_m=along_leaf,
            component_norm_share_m=components@(component/np.linalg.norm(component))
                if np.linalg.norm(component) else np.zeros(len(block)),
            summed_component_m=components.sum(0),summed_mid_action_m=mid.sum(0),radii_m=radii)
        for observed, expected in [(sum(scalar),remainder),(sum(along_leaf),along),
                                   (components.sum(0),component),(mid.sum(0),action)]:
            np.testing.assert_allclose(observed,expected,rtol=1e-11,atol=1e-21)
    return action,""")
    scope = dict(vars(previous.p16))
    exec(compile(source, str(HERE/"face.py"), "exec"), scope)
    return scope["paired"](red, refined, displacements, radii, leaf_box, kt, km, A, ix, direction, dest)


def run(sign, destination, tight):
    def nominal(*args):
        step = int(args[-1].stem.split("-")[1])
        # Use the whole ranked cohort at 832. At 833 use the same 328 leaves
        # as p18 so that the measured change comes only from the previous radius.
        previous.negative_group = sys.maxsize if tight and step == 832 else 328
        return previous.integrate(*args)
    previous.p15.integrate = nominal
    previous.p16.paired = measured
    previous.p16.run(sign, destination, "projected")


if __name__ == "__main__":
    run(sys.argv[1], Path(sys.argv[2]), "--tight" in sys.argv[3:])
