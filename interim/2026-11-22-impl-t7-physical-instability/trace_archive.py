"""Losslessly transcode generated gzip JSON for the repository evidence archive."""

import argparse
import gzip
import hashlib
import json
import lzma
from pathlib import Path
import shutil


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def transcode(src, dst, opening, writing):
    with opening(src) as inp, writing(dst) as out:
        shutil.copyfileobj(inp, out, 1024**2)


def compress(source, destination):
    destination.mkdir(parents=True, exist_ok=True)
    entries = {}
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(source)
        out = destination / relative
        out.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            packed = out.with_suffix(".xz")
            if not packed.exists():
                transcode(
                    path, packed, gzip.open, lambda p: lzma.open(p, "xb", preset=1)
                )
            entries[relative.as_posix()] = dict(
                gzip_sha256=sha(path), xz_sha256=sha(packed)
            )
        else:
            shutil.copyfile(path, out)
    (destination / "trace-encoding.json").write_text(
        json.dumps(entries, indent=2) + "\n", encoding="utf-8"
    )


def restore(folder):
    entries = json.loads((folder / "trace-encoding.json").read_bytes())
    for name, has in entries.items():
        out = folder / name
        packed = out.with_suffix(".xz")
        if sha(packed) != has["xz_sha256"]:
            raise ValueError(f"XZ digest differs: {packed}")
        if out.exists():
            if sha(out) != has["gzip_sha256"]:
                raise ValueError(f"Existing trace differs: {out}")
            continue
        with lzma.open(packed) as src, out.open("xb") as stream:
            with gzip.GzipFile(
                filename="", fileobj=stream, mode="wb", mtime=0, compresslevel=1
            ) as dst:
                shutil.copyfileobj(src, dst, 1024**2)
        if sha(out) != has["gzip_sha256"]:
            raise ValueError(f"Gzip implementation or trace differs: {out}")
    return len(entries)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--compress", nargs=2, type=Path, metavar=("SOURCE", "DEST"))
    g.add_argument("--restore", type=Path)
    a = p.parse_args()
    if a.compress:
        compress(*a.compress)
    else:
        print(restore(a.restore))
