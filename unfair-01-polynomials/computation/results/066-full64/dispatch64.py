"""Run exact degree-64 shards and require complete, disjoint coverage."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import time

MASK_COUNT = 1 << 23
SHARD_SIZE = 1 << 16


def run_shard(binary, directory, start):
    stop = start + SHARD_SIZE
    began = time.monotonic()
    with subprocess.Popen(
        [str(binary), str(start), str(stop)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ) as proc:
        try:
            output, errors = proc.communicate(timeout=180)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.communicate()
            raise
        assert proc.returncode == 0 and not errors, (start, errors)
    row = json.loads(output)
    assert (row["start"], row["stop"]) == (start, stop)
    assert 0 <= len(row["hits"]) <= row["first_prime_pass"] <= SHARD_SIZE
    for hit in row["hits"]:
        assert start <= hit["mask"] < stop
        assert min(hit["gcd163"], hit["gcd1279"]) >= 7
    row["seconds"] = time.monotonic() - began
    row["pid"] = proc.pid
    checkpoint = directory / f"{start:08d}.json"
    with checkpoint.open("x") as handle:
        json.dump(row, handle, sort_keys=True)
        handle.write("\n")
    return row


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    out = args.output
    out.mkdir(exist_ok=False, parents=True)
    work = range(0, MASK_COUNT, SHARD_SIZE)
    began = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(lambda start: run_shard(args.binary.resolve(), out, start), work))
    assert [row["start"] for row in rows] == list(work)
    report = {
        "status": "complete",
        "masks": MASK_COUNT,
        "shard_size": SHARD_SIZE,
        "workers": args.workers,
        "shards": len(rows),
        "first_prime_pass": sum(row["first_prime_pass"] for row in rows),
        "hits": [hit for row in rows for hit in row["hits"]],
        "seconds": time.monotonic() - began,
    }
    with (out / "summary.json").open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
