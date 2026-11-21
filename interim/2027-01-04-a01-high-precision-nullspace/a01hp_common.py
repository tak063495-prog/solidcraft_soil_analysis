"""Standard-library-only evidence I/O; no high-precision dependency in CI."""
import json
from pathlib import Path

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(value,f,indent=2,ensure_ascii=False);f.write("\n")

def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
