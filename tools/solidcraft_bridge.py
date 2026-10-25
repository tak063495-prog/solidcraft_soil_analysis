"""Snapshot the running SolidCraft design over its authenticated MCP bridge."""
import argparse
import json
from pathlib import Path
import queue
import subprocess
import threading


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, required=True)
    parser.add_argument("--connect", default="127.0.0.1:7878")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    step = output / "model.step"
    if step.exists():
        parser.error("model.step already exists")
    # stderr goes to a file, so a full error pipe cannot deadlock a long export.
    with (output / "solidcraft-mcp.stderr.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen([str(args.cli), "mcp", "--connect", args.connect,
                                    "--allow-path", str(output)], stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=log, text=True, encoding="utf-8")
        replies = queue.Queue()

        def reader():
            for line in process.stdout:
                replies.put(line)
            replies.put(None)

        threading.Thread(target=reader, daemon=True).start()

        def send(payload):
            process.stdin.write(json.dumps(payload) + "\n")
            process.stdin.flush()

        def request(identifier, method, params):
            send({"jsonrpc": "2.0", "id": identifier, "method": method, "params": params})
            while True:
                line = replies.get(timeout=120)
                if line is None:
                    raise RuntimeError("SolidCraft MCP exited: see solidcraft-mcp.stderr.log")
                response = json.loads(line)
                if response.get("id") != identifier:
                    continue
                result = response.get("result", {})
                if "error" in response or result.get("isError"):
                    raise RuntimeError(json.dumps(response, ensure_ascii=False))
                return result

        try:
            request(1, "initialize", {"protocolVersion": "2025-06-18", "capabilities": {},
                                       "clientInfo": {"name": "soil-fem-3d", "version": "0.1.0"}})
            send({"jsonrpc": "2.0", "method": "notifications/initialized"})
            before = request(2, "tools/call", {"name": "inspect_design", "arguments": {}})
            export = request(3, "tools/call", {"name": "export", "arguments": {"path": str(step), "format": "step"}})
            after = request(4, "tools/call", {"name": "inspect_design", "arguments": {}})
            if not step.is_file() or step.stat().st_size == 0:
                raise RuntimeError("STEP export missing")
            # Detect edits during capture; mesh only a stable model snapshot.
            if before != after:
                raise RuntimeError("design changed during export; retry in a new run directory")
            (output / "solidcraft-inspect.json").write_text(json.dumps(before, indent=2), encoding="utf-8")
            (output / "solidcraft-export.json").write_text(json.dumps(export, indent=2), encoding="utf-8")
            print(step)
        finally:
            process.stdin.close()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=5)


if __name__ == "__main__":
    main()
