#!/usr/bin/env python3
"""Sequential exact n=7 scan of D4-inequivalent four-point critical clusters."""
import json
from pathlib import Path
import subprocess
import sys

n = 7
base = [(1, 0), (2, 0), (3, 0), (1, 1)]


def mask(cells):
    return sum(1 << (y * n + x) for x, y in cells)


def transform(cells, op):
    out = []
    for x, y in cells:
        if op >= 4:
            x = n - 1 - x
        for _ in range(op % 4):
            x, y = n - 1 - y, x
        out.append((x, y))
    return tuple(sorted(out))


all_masks = set()
for op in range(8):
    shape = transform(base, op)
    minx, maxx = min(x for x, _ in shape), max(x for x, _ in shape)
    miny, maxy = min(y for _, y in shape), max(y for _, y in shape)
    for dx in range(-minx, n - maxx):
        for dy in range(-miny, n - maxy):
            all_masks.add(mask([(x + dx, y + dy) for x, y in shape]))


def d4mask(m, op):
    return mask(transform([(v % n, v // n) for v in range(n * n) if m >> v & 1], op))


reps = sorted({min(d4mask(m, op) for op in range(8)) for m in all_masks})


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: run_n7_critical_orbits.py BINARY OUTPUT_DIR")
    binary, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    if out.exists():
        raise SystemExit(f"refusing existing {out}")
    out.mkdir(parents=True)
    manifest = {"n": n, "shape": base, "all_transforms_and_translations": len(all_masks),
                "d4_orbit_representatives": [hex(x) for x in reps], "runs": []}
    for i, m in enumerate(reps):
        stem = f"orbit-{i:02d}-{m:x}"
        command = ["/usr/bin/time", "-l", "-o", str(out / f"{stem}.time"), str(binary), "7", hex(m)]
        with (out / f"{stem}.json").open("wb") as fo, (out / f"{stem}.stderr").open("wb") as fe:
            p = subprocess.run(command, stdout=fo, stderr=fe, check=False)
        item = {"index": i, "mask": hex(m), "command": command, "exit_code": p.returncode,
                "json": f"{stem}.json", "time": f"{stem}.time", "stderr": f"{stem}.stderr"}
        if p.returncode == 0:
            try:
                x = json.loads((out / item["json"]).read_text())
                item["candidate_holds"] = x["candidate_holds"]
                item["independent_sets"] = x["independent_sets"]
            except Exception as exc:
                item["parse_error"] = repr(exc)
                item["exit_code"] = 99
        manifest["runs"].append(item)
        (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if item["exit_code"]:
            raise SystemExit(item["exit_code"])
    manifest["status"] = "completed"
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
