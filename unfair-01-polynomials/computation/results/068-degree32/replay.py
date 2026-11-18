"""Restore and check the exact degree-32 witness slices."""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "degree32-archive"
SOURCES = [HERE / "source", HERE / "round065-gap7-next", HERE / "round066-gap7-next"]


def digest(data):
    return sha256(data).hexdigest()


def restore(parts, destination, total_sha):
    total = sha256()
    with destination.open("xb") as output:
        for part in parts:
            compressed = (ARCHIVE / part["file"]).read_bytes()
            assert len(compressed) == part["compressed_bytes"]
            assert digest(compressed) == part["compressed_sha256"]
            raw = gzip.decompress(compressed)
            assert len(raw) == part["raw_bytes"] and digest(raw) == part["raw_sha256"]
            output.write(raw)
            total.update(raw)
    assert total.hexdigest() == total_sha


def run(out, label, module, *args):
    command = [sys.executable, "-B", str(SOURCES[0] / (module+".py")), *map(str, args)]
    with (out / (label+".log")).open("xb") as log:
        subprocess.run(command, check=True, cwd=out, stdout=log, stderr=subprocess.STDOUT)


def check_part(out, census, i, case):
    raw = out / f"raw-{i}.jsonl"
    restore(case["parts"], raw, case["raw_sha256"])
    replay = out / f"verified-{i}.jsonl"
    run(out, f"slice-{i}", "certify32", "--census", census, "--start", case["start"],
        "--stop", case["stop"], "--verify", raw, "--out", replay)
    with replay.open() as stream:
        json.loads(next(stream))
        count = 0
        for line in stream:
            final = json.loads(line)
            count += 1
    assert final["status"] == "complete"
    assert (final["start"], final["stop"]) == (case["start"], case["stop"])
    assert count-1 == final["stop"]-final["start"]
    assert final["counts"] == case["counts"]
    return final


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    m = json.loads((ARCHIVE / "verification.json").read_bytes())
    for name, expected in m["sources_sha256"].items():
        files = [p/name for p in SOURCES if (p/name).is_file()]
        assert len(files) == 1
        assert digest(files[0].read_bytes()) == expected
    census = out / "census.json"
    restore(m["census"]["parts"], census, m["census"]["sha256"])
    audit = out / "census-audit.json"
    real = out / "real-check.json"
    cubic = out / "cubic-check.json"
    obstruction = HERE / "degree32-first-unresolved.json"
    ordered=sorted((part["start"],part["stop"]) for part in m["slices"])
    assert ordered[0][0]==0 and ordered[-1][-1]==m["counts"]["shapes"]["candidates"]
    for left,right in zip(ordered,ordered[1:]):
        assert left[0]<left[1]==right[0]
    assert ordered[-1][0]<ordered[-1][1]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        fa = pool.submit(run, out, "census", "audit_census32", "--census", census, "--out", audit)
        fr = pool.submit(run, out, "real", "verify_real32", "--obstruction", obstruction, "--out", real)
        fc = pool.submit(run, out, "cubic", "verify_cubic32", "--census", census, "--out", cubic)
        slices = [pool.submit(check_part, out, census, i, case) for i, case in enumerate(m["slices"])]
        fa.result()
        fr.result()
        fc.result()
        finished = [piece.result() for piece in slices]
    counts = Counter()
    for case in finished:
        counts.update(case["counts"])
    assert counts == m["counts"]["modular"]
    assert digest(real.read_bytes()) == m["real_check"]["sha256"]
    assert digest(cubic.read_bytes()) == m["cubic_check"]["sha256"]
    report = dict(status="verified", workers=args.workers, counts=dict(counts), slices=finished,
                  resolved_in_real=m["real_check"]["resolved"])
    with (out / "receipt.json").open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(report["status"], counts["unit"], len(report["resolved_in_real"]))


if __name__ == "__main__":
    main()
