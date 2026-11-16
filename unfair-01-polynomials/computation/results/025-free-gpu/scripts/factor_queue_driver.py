#!/usr/bin/env python3
"""Prepare or run bounded, recoverable MPS factor-search batches.

Dry-run is the default.  --execute and exactly one --job are required to
launch a child.  The frozen worker writes frequent numerical checkpoints;
this driver supplies a four-hour wall-time boundary and records PID/release
evidence.  It never claims an optimisation failure is a proof.
"""
from __future__ import annotations

import argparse, hashlib, json, os, signal, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKER = HERE / "gpu_factor_search_frozen.py"
JOBS = (
    {"name": "m20n80d1", "m": 20, "n": 80, "fractional_index": 1, "starts": 16384, "margin": 0.05, "seed": 2026111601},
    {"name": "m24n64d2", "m": 24, "n": 64, "fractional_index": 2, "starts": 16384, "margin": 0.05, "seed": 2026111602},
    {"name": "m28n60d3", "m": 28, "n": 60, "fractional_index": 3, "starts": 16384, "margin": 0.05, "seed": 2026111603},
)

def digest(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def atomic(path: Path, value: object) -> None:
    temp = path.with_suffix(path.suffix + ".tmp"); temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n"); temp.replace(path)
def command(job: dict, run_dir: Path, steps: int) -> list[str]:
    return [sys.executable, str(WORKER), "--m", str(job["m"]), "--n", str(job["n"]), "--fractional-index", str(job["fractional_index"]), "--starts", str(job["starts"]), "--margin", str(job["margin"]), "--seed", str(job["seed"]), "--steps", str(steps), "--restart-every", "2000", "--out", str(run_dir)]
def stop_owned_child(proc: subprocess.Popen, grace_seconds: int) -> str:
    os.killpg(proc.pid, signal.SIGINT)
    try: proc.wait(timeout=grace_seconds); return "SIGINT"
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGTERM)
        try: proc.wait(timeout=grace_seconds); return "SIGTERM"
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL); proc.wait(); return "SIGKILL"

def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--run-root", type=Path, default=HERE / "runs")
    ap.add_argument("--steps", type=int, default=500000); ap.add_argument("--time-limit-seconds", type=int, default=4 * 3600)
    ap.add_argument("--poll-seconds", type=int, default=30); ap.add_argument("--job", choices=[j["name"] for j in JOBS])
    ap.add_argument("--execute", action="store_true", help="launch only the named owned child")
    args = ap.parse_args()
    if not Path("/Volumes/data").is_mount() or not WORKER.is_file() or min(args.steps, args.time_limit_seconds, args.poll_seconds) <= 0: raise SystemExit("invalid mount, worker, or bounds")
    if args.execute != bool(args.job): raise SystemExit("--execute requires exactly one --job; dry-run launches nothing")
    args.run_root.mkdir(parents=True, exist_ok=True); queue_path = args.run_root / "queue.json"
    scope = ("Numerical sparse-low-prefix subclass only: both factors force coefficients 1..d-1 to zero, "
             "coefficient d is fractional in [margin,1-margin], and every later coefficient is independently optimized then projected only to [0,1]. "
             "Thus gaps after d are allowed, but this is not the unrestricted first-fractional-index model or a proof-level exclusion.")
    queue = {"format":"bounded-factor-search-queue-v1", "worker":str(WORKER), "worker_sha256":digest(WORKER), "scope":scope,
             "time_limit_seconds_per_child":args.time_limit_seconds, "jobs":[{**j,"command":command(j,args.run_root/j["name"],args.steps)} for j in JOBS], "status":"prepared"}
    atomic(queue_path,queue)
    if not args.execute:
        print(json.dumps({"status":"prepared","queue":str(queue_path),"jobs":[j["name"] for j in JOBS]})); return
    job = next(j for j in JOBS if j["name"] == args.job); run_dir=args.run_root/job["name"]
    if run_dir.exists() and any(run_dir.iterdir()): raise SystemExit(f"refusing to overwrite nonempty run directory {run_dir}")
    run_dir.mkdir(exist_ok=True); log=(run_dir/"driver-child.log").open("w")
    proc=subprocess.Popen(command(job,run_dir,args.steps),stdout=log,stderr=subprocess.STDOUT,start_new_session=True); started=time.monotonic()
    entry=next(j for j in queue["jobs"] if j["name"]==job["name"]); entry.update({"pid":proc.pid,"status":"running"}); queue["status"]="running"; atomic(queue_path,queue)
    try:
        while proc.poll() is None and time.monotonic()-started < args.time_limit_seconds: time.sleep(args.poll_seconds)
        state="exited" if proc.poll() is not None else "time_boundary_"+stop_owned_child(proc,120)
    except BaseException:
        if proc.poll() is None: stop_owned_child(proc,120)
        entry.update({"status":"driver_interrupted_after_owned_child_recovery","returncode":proc.returncode}); atomic(queue_path,queue); raise
    finally: log.close()
    entry.update({"status":state,"returncode":proc.returncode,"elapsed_seconds":time.monotonic()-started,"checkpoint_exists":(run_dir/"checkpoint.json").is_file(),"release_exists":(run_dir/"release.json").is_file()}); queue["status"]="child_finished"; atomic(queue_path,queue)

if __name__ == "__main__": main()
