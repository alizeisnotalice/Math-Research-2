#!/usr/bin/env python3
"""Float32 CPU/MPS calibration only; C++ is the exact certification path."""
import argparse
import json
import os
import platform
import time
os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "0"
import torch

parser = argparse.ArgumentParser()
parser.add_argument("--batch", type=int, default=131072)
parser.add_argument("--reps", type=int, default=5)
parser.add_argument("--out", required=True)
args = parser.parse_args()
torch.manual_seed(20261116)
choices = torch.randint(0, 8, (args.batch, 8), dtype=torch.int64)
choices[:, 0] %= 4
phase = torch.tensor([[1 if (a & u).bit_count() % 2 == 0 else -1 for u in range(8)] for a in range(8)], dtype=torch.float32)

def run(device: str):
    z = choices.to(device)
    b = z & 3
    signs = torch.where(z < 4, 1., -1.)
    signs[:, 0] = 1.
    q = phase.to(device)
    values = [(signs * (b == c).float()) @ q.T for c in range(4)]
    return torch.stack(values, dim=2)

reference = run("cpu")
report = {"seed":20261116,"precision":"float32", "comparison":"all Walsh entries, CPU vs device", "batch": args.batch, "reps": args.reps, "torch": torch.__version__, "platform": platform.platform(), "runs": {}}
devices = [("cpu", "cpu")]
if torch.backends.mps.is_available():
    devices.append(("mps", "mps"))
for name, device in devices:
    run(device)
    if device == "mps":
        torch.mps.synchronize()
    times, values = [], []
    for _ in range(args.reps):
        began = time.perf_counter()
        value = run(device)
        if device == "mps":
            torch.mps.synchronize()
        host_value = value.cpu()
        times.append(time.perf_counter() - began)
        if not torch.equal(host_value, reference):
            raise RuntimeError("entrywise CPU/device mismatch")
        values.append(int(host_value.ne(0).sum()))
    report["runs"][name] = {"device": device, "seconds": times, "median_seconds": sorted(times)[len(times)//2], "values": values}
if "mps" in report["runs"] and set(report["runs"]["cpu"]["values"]) != set(report["runs"]["mps"]["values"]):
    raise RuntimeError("CPU/MPS mismatch")
open(args.out, "w").write(json.dumps(report, indent=2) + "\n")
