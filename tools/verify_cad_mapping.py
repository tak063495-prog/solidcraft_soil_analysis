"""Exercise named-body export failure cases and fragment-ancestry ordering."""
import argparse
import copy
import json
from pathlib import Path
from cad_layers import snapshot
from mesh_layers import generate


def verify(cli, base, output):
    output.mkdir(parents=True,exist_ok=False)
    report = {}
    ground = json.loads((base/"ground.json").read_text(encoding="utf-8"))
    invalid = copy.deepcopy(ground)
    invalid["layers"][0]["body"] = "UnknownBody"
    path = output/"invalid-ground.json"
    path.write_text(json.dumps(invalid),encoding="utf-8")
    try:
        snapshot(cli,output/"unknown-body",path,Path("examples/layers.json"))
        raise AssertionError("unknown body was accepted")
    except ValueError as error:
        assert "mapping mismatch" in str(error),error
        report["unknown_body_mapping"] = "rejected"
    assert not (output/"unknown-body/bodies.json").exists()

    model = json.loads(Path("examples/layers.json").read_text(encoding="utf-8"))
    model["commands"][1]["params"]["corner"][2] = 4000
    path = output/"overlap.json"
    path.write_text(json.dumps(model),encoding="utf-8")
    snapshot(cli,output/"overlap",base/"ground.json",path)
    try:
        generate(output/"overlap/bodies.json",output/"overlap/mesh.json",2.)
        raise AssertionError("overlapping solids were accepted")
    except ValueError as error:
        assert "overlap in positive volume" in str(error),error
        report["overlapping_bodies"] = "rejected"
    assert not (output/"overlap/mesh.json").exists()

    manifest = json.loads((base/"bodies.json").read_text(encoding="utf-8"))
    for r in manifest["regions"]:
        r["step"] = str((base/r["step"]).resolve())
    manifest["regions"].reverse()
    path = output/"reversed-bodies.json"
    path.write_text(json.dumps(manifest),encoding="utf-8")
    generate(path,output/"reversed-mesh.json",2.)
    result = json.loads((output/"reversed-mesh.json").read_text(encoding="utf-8"))
    for i,tet in enumerate(result["tetrahedra"]):
        z = sum(result["nodes"][n][2] for n in tet)/4
        assert result["material_ids"][i] == (1 if z<5 else 2)
    report["reversed_step_import_order"] = "material identity preserved"

    manifest["regions"][0]["source_sha256"] = "0"*64
    path = output/"corrupt-manifest.json"
    path.write_text(json.dumps(manifest),encoding="utf-8")
    try:
        generate(path,output/"corrupt-mesh.json",2.)
        raise AssertionError("STEP hash mismatch was accepted")
    except ValueError as error:
        assert "hash mismatch" in str(error),error
        report["step_hash_mismatch"] = "rejected"
    (output/"cad-mapping-verification.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli",type=Path,required=True)
    parser.add_argument("--base-run",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    a = parser.parse_args()
    verify(a.cli.resolve(),a.base_run.resolve(),a.output_dir.resolve())
