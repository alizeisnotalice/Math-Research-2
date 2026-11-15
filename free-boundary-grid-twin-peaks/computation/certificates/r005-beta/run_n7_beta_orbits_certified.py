#!/usr/bin/env python3
"""Run the 12 n=7 D4 edge orbits sequentially with complete per-run records.

Each subprocess itself uses the C++ program's 12 workers, but this driver
waits for it to exit before starting the next orbit.  It creates a fresh
output directory and refuses to overwrite it.
"""
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--orbits", type=int, default=12)
    args = ap.parse_args()
    binary, source, root = args.binary.resolve(), args.source.resolve(), args.output_dir.resolve()
    if not binary.is_file() or not source.is_file():
        raise SystemExit("binary or source is missing")
    if root.exists():
        raise SystemExit(f"refusing to overwrite existing {root}")
    root.mkdir(parents=True)
    manifest = {"n": 7, "orbits_requested": args.orbits,
                "source": str(source), "source_sha256_before": sha256(source),
                "binary": str(binary), "binary_sha256_before": sha256(binary),
                "runs": []}
    for orbit in range(args.orbits):
        stem = f"orbit-{orbit:02d}"
        out, err, timing = (root / f"{stem}.json", root / f"{stem}.stderr", root / f"{stem}.time")
        command = ["/usr/bin/time", "-l", "-o", str(timing), str(binary), "7", "--orbit", str(orbit)]
        start = now()
        with out.open("wb") as f_out, err.open("wb") as f_err:
            proc = subprocess.run(command, stdout=f_out, stderr=f_err, check=False)
        end = now()
        item = {"orbit": orbit, "command": command, "started_at": start,
                "finished_at": end, "exit_code": proc.returncode,
                "json": out.name, "stderr": err.name, "time": timing.name}
        if proc.returncode == 0:
            try:
                result = json.loads(out.read_text())
                ok = (result.get("n") == 7 and result.get("selected_orbit") == orbit
                      and len(result.get("edge_orbits", [])) == 1
                      and result.get("directed_edge_coverage")
                         == 2 * result["edge_orbits"][0]["multiplicity"])
                item["post_exit_json_check"] = bool(ok)
                if not ok:
                    item["exit_code"] = 98
            except Exception as exc:  # preserve the process files, then stop
                item["post_exit_json_check"] = False
                item["json_error"] = repr(exc)
                item["exit_code"] = 99
        manifest["runs"].append(item)
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if item["exit_code"] != 0:
            raise SystemExit(item["exit_code"])
    manifest["source_sha256_after"] = sha256(source)
    manifest["binary_sha256_after"] = sha256(binary)
    manifest["status"] = "all_orbits_completed"
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
