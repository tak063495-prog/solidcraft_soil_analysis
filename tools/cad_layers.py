"""Body-resolved SolidCraft snapshot: STEP ordering never determines material identity."""
import argparse
import hashlib
import json
from pathlib import Path
import queue
import subprocess
import threading


def step_volume(path):
    """SolidCraft measures by fine tessellation; STEP is the analysis solid."""
    import gmsh
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal',0)
        gmsh.option.setString('Geometry.OCCTargetUnit','MM')
        volumes=[x for x in gmsh.model.occ.importShapes(str(path.resolve())) if x[0]==3]
        if len(volumes)!=1:raise ValueError('each named body must export one closed STEP solid')
        return gmsh.model.occ.getMass(*volumes[0])*1e-9
    finally:gmsh.finalize()


def snapshot(cli, output, ground, design=None, connect=None):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    target = output / "bodies.json"
    if target.exists():
        raise ValueError("body snapshot already exists")
    case_bytes = ground.read_bytes()
    config = json.loads(case_bytes)
    if config.get("schema") != "soil-fem-ground/1" or config.get("units") != "m-kN-kPa":
        raise ValueError("ground schema/units mismatch")
    bindings=config.get('bodies',[]) if config.get('initialization',{}).get('method')=='equilibrate' else config.get('layers',[])
    mapping = {x["body"]: x["material_id"] for x in bindings}
    tools = config.get("excavation_bodies", [])
    if len(set(tools)) != len(tools) or set(tools) & set(mapping):
        raise ValueError("excavation tools must have unique names distinct from material bodies")
    if not mapping or len(mapping) != len(bindings):
        raise ValueError("duplicate or missing body mapping")
    command = [str(cli), "mcp", "--allow-path", str(output)]
    if connect:
        command += ["--connect", connect]
    else:
        command += ["--in", str(design.resolve()), "--allow-path", str(design.resolve().parent)]
    with (output / "cad.stderr.log").open("x", encoding="utf-8") as log:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=log, text=True, encoding="utf-8")
        lines = queue.Queue()

        def read():
            for line in process.stdout:
                lines.put(line)
            lines.put(None)

        threading.Thread(target=read, daemon=True).start()

        def send(message):
            process.stdin.write(json.dumps(message) + "\n")
            process.stdin.flush()

        def request(i, method, params):
            send({"jsonrpc":"2.0", "id":i, "method":method, "params":params})
            while True:
                line = lines.get(timeout=120)
                if line is None:
                    raise RuntimeError("SolidCraft exited; inspect cad.stderr.log")
                reply = json.loads(line)
                if reply.get("id") != i:
                    continue
                result = reply.get("result", {})
                if "error" in reply or result.get("isError"):
                    raise RuntimeError(json.dumps(reply))
                return result

        def call(i, name, arguments):
            result = request(i, "tools/call", {"name":name, "arguments":arguments})
            if "structuredContent" in result:
                return result["structuredContent"]
            return json.loads(next(x["text"] for x in result["content"] if x["type"] == "text"))

        try:
            request(1,"initialize",{"protocolVersion":"2025-06-18","capabilities":{},
                                    "clientInfo":{"name":"soil-fem-layers","version":"0.2.0"}})
            send({"jsonrpc":"2.0","method":"notifications/initialized"})
            before = call(2,"inspect_design",{})
            bodies = before["bodies"]
            names = [b["name"] for b in bodies]
            if len(names) != len(set(names)) or set(names) != set(mapping) | set(tools):
                raise ValueError(f"CAD/body mapping mismatch: CAD={names}, configured={list(mapping)}")
            if any(f.get("error") for f in before.get("timeline", [])):
                raise ValueError("CAD timeline has errors")
            if before.get("units") != "mm":
                raise ValueError("this adapter requires a millimetre CAD design")
            records, tool_records = [], []
            for i, body in enumerate(sorted(bodies, key=lambda x: x["name"])):
                path = output / f"body-{i+1:04}.step"
                if path.exists():
                    raise ValueError(f"snapshot output exists: {path}")
                call(10+i,"export",{"path":str(path),"format":"step","bodies":[body["name"]]})
                measured=body["volume_mm3"]*1e-9
                exact=step_volume(path)
                relative_difference=abs(exact-measured)/max(abs(exact),abs(measured),1e-30)
                if exact<=0 or relative_difference>0.002:
                    raise ValueError(f'STEP/display volume mismatch for {body["name"]}: reported={measured}, STEP={exact}')
                record = {"id":i+1,"body":body["name"],
                                "step":path.name,"source_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
                                "bounds_m":[[v*0.001 for v in body["bbox"][key]] for key in ("min","max")],
                                "volume_m3":exact,"reported_tessellated_volume_m3":measured,
                                "step_vs_reported_relative_difference":relative_difference}
                if body["name"] in mapping:
                    record["material_id"] = mapping[body["name"]]
                    records.append(record)
                else:
                    tool_records.append(record)
            after = call(100000,"inspect_design",{})
            if before != after:
                raise ValueError("CAD changed while exporting bodies; retry in a fresh directory")
            if ground.read_bytes() != case_bytes:
                raise ValueError("ground configuration changed during capture")
            manifest = {"schema":"soil-fem-cad-bodies/1","units":"m-kN-kPa",
                        "ground_sha256":hashlib.sha256(case_bytes).hexdigest(),"regions":records}
            if tool_records:
                manifest["excavation_tools"] = tool_records
            manifest["volume_measurement"]={"analysis":"STEP B-rep","reported":"SolidCraft fine tessellation","maximum_relative_difference":0.002}
            target.write_text(json.dumps(manifest,indent=2,allow_nan=False),encoding="utf-8")
            (output / "cad-inspect.json").write_text(json.dumps(before,indent=2),encoding="utf-8")
            print(f"SolidCraft: captured {len(records)} material bodies and {len(tool_records)} excavation tools")
        finally:
            process.stdin.close()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--ground",type=Path,required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--design",type=Path)
    source.add_argument("--connect")
    a = parser.parse_args()
    snapshot(a.cli,a.output_dir,a.ground,a.design,a.connect)
