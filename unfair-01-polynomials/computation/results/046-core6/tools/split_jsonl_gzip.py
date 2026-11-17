#!/usr/bin/env python3
"""Split a gzip JSONL audit trail on line boundaries, retaining provenance."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("out_dir", type=Path)
    p.add_argument("--stem", default="single-terminal-witnesses")
    p.add_argument("--plain-part-bytes", type=int, default=8 * 1024 * 1024)
    args = p.parse_args()
    if args.plain_part_bytes <= 0:
        raise SystemExit("--plain-part-bytes must be positive")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_sha = digest(args.source)
    source_plain_sha = hashlib.sha256()
    parts: list[dict[str, object]] = []
    line_count = 0
    part_index = 0
    plain_in_part = 0
    out_f = None
    out_gz = None
    part_path: Path | None = None

    def close_part() -> None:
        nonlocal out_f, out_gz, part_path, plain_in_part
        if out_gz is None:
            return
        out_gz.close()
        out_f.close()
        assert part_path is not None
        size = part_path.stat().st_size
        if size >= 10 * 1024 * 1024:
            raise RuntimeError(f"part exceeds 10 MiB: {part_path} ({size} bytes)")
        parts.append({
            "name": part_path.name,
            "compressed_bytes": size,
            "compressed_sha256": digest(part_path),
            "uncompressed_bytes": plain_in_part,
        })
        out_f = out_gz = part_path = None
        plain_in_part = 0

    with gzip.open(args.source, "rb") as inp:
        for line in inp:
            if out_gz is None or (plain_in_part and plain_in_part + len(line) > args.plain_part_bytes):
                close_part()
                part_path = args.out_dir / f"{args.stem}.part-{part_index:03d}.jsonl.gz"
                part_index += 1
                out_f = part_path.open("wb")
                out_gz = gzip.GzipFile(filename="", mode="wb", fileobj=out_f, mtime=0)
            out_gz.write(line)
            source_plain_sha.update(line)
            plain_in_part += len(line)
            line_count += 1
    close_part()
    meta = {
        "format": "line-preserving gzip split; concatenate decompressed parts in listed order",
        "source_raw_filename": args.source.name,
        "source_raw_sha256": raw_sha,
        "source_uncompressed_sha256": source_plain_sha.hexdigest(),
        "line_count": line_count,
        "parts": parts,
    }
    (args.out_dir / f"{args.stem}-split.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({"parts": len(parts), "line_count": line_count, "raw_sha256": raw_sha}, indent=2))


if __name__ == "__main__":
    main()
