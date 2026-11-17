#!/usr/bin/env python3
"""Local, low-memory orchestrator for the authorized N=9 full scan.

It creates one fresh run directory, compiles with the shared module cache,
runs both gates and a counted benchmark prefix, then keeps up to twelve C++
workers busy with dynamically issued small chunks.  It never touches another
process or an existing run directory.
"""
from __future__ import annotations

from collections import deque
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, os, subprocess, sys, time

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "bipartite_prefix_dp.cpp"
# Portable public-package adaptation: callers supply a writable shared cache.
CACHE = Path(os.environ["CLANG_MODULE_CACHE_PATH"])
TOTALS = {1: 256, 2: 16384, 3: 262144, 4: 1048576}
WORKERS, CHUNK = 12, 4096


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic(path: Path, obj: object) -> None:
    if path.exists():
        raise FileExistsError(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def invoke(command: list[str], cwd: Path, stdout: Path, stderr: Path, env: dict[str, str]) -> dict:
    started, t0 = utc(), time.perf_counter()
    with stdout.open("wb") as out, stderr.open("wb") as err:
        proc = subprocess.run(command, cwd=cwd, stdout=out, stderr=err, env=env, check=False)
    return {"command": command, "exit_code": proc.returncode, "utc_started": started,
            "utc_finished": utc(), "wall_seconds": time.perf_counter()-t0,
            "stdout_sha256": sha(stdout), "stderr_sha256": sha(stderr)}


def task_list() -> deque[dict]:
    tasks: deque[dict] = deque()
    # The benchmark below owns [0,1024) for left=4 and is counted exactly once.
    for left, start in ((1, 0), (2, 0), (3, 0), (4, 1024)):
        total = TOTALS[left]
        while start < total:
            stop = min(start + CHUNK, total)
            tasks.append({"left": left, "start": start, "stop": stop,
                          "id": f"prod-l{left}-{start:010d}-{stop:010d}"})
            start = stop
    return tasks


def checkpoint_name(task: dict) -> str:
    return f"shard-l{task['left']}-{task['id']}.json"


def direct_verify(binary: Path, run: Path, source_sha: str, task: dict, checkpoint: Path, env: dict[str, str]) -> dict:
    row = json.loads(checkpoint.read_text())
    assert row["status"] == "counterexample"
    edge_mask = str(row["counterexample"]["edge_mask"])
    logbase = run / "counterexample" / task["id"]
    logbase.parent.mkdir(exist_ok=True)
    cmd = ["nice", "-n", "10", str(binary), "--mode", "verify-edge-mask-r27", "--left-size", str(task["left"]),
           "--edge-mask", edge_mask, "--out-dir", str(logbase.parent), "--source-sha", source_sha]
    receipt = invoke(cmd, run, logbase.with_suffix(".stdout"), logbase.with_suffix(".stderr"), env)
    receipt.update({"task": task, "checkpoint": str(checkpoint), "edge_mask": edge_mask})
    atomic(logbase.with_suffix(".receipt.json"), receipt)
    return receipt


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    run = ROOT / f"local-cpu-r034-{stamp}"
    if run.exists():
        raise SystemExit(f"refusing existing run directory: {run}")
    run.mkdir(); (run / "results").mkdir(); (run / "gates").mkdir(); (run / "logs").mkdir()
    CACHE.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ); env["CLANG_MODULE_CACHE_PATH"] = str(CACHE)
    binary = run / "bipartite_prefix_dp"
    source_sha = sha(SOURCE)
    atomic(run / "started.json", {"utc_started": utc(), "source": str(SOURCE), "source_sha256": source_sha,
                                   "cache": str(CACHE), "workers": WORKERS, "chunk": CHUNK,
                                   "benchmark_counted_range": {"left": 4, "start": 0, "stop": 1024}})
    compile_cmd = ["c++", "-O3", "-std=c++20", "-Wall", "-Wextra", "-Werror", "-pedantic", str(SOURCE), "-o", str(binary)]
    compile_receipt = invoke(compile_cmd, run, run / "compile.stdout", run / "compile.stderr", env)
    compile_receipt["binary_sha256"] = sha(binary) if binary.exists() else None
    atomic(run / "compile.receipt.json", compile_receipt)
    if compile_receipt["exit_code"]:
        raise SystemExit("compile failed; see receipt")
    for mode in ("gate-n7-r45", "gate-n9-direct"):
        receipt = invoke([str(binary), "--mode", mode, "--out-dir", str(run / "gates"), "--source-sha", source_sha],
                         run, run / "gates" / f"{mode}.stdout", run / "gates" / f"{mode}.stderr", env)
        atomic(run / "gates" / f"{mode}.receipt.json", receipt)
        if receipt["exit_code"]:
            raise SystemExit(f"{mode} failed or found a counterexample; production not launched")
    # The benchmark is a real production coverage prefix, not disposable work.
    bench = {"left": 4, "start": 0, "stop": 1024, "id": "benchmark-l4-0000000000-0000001024"}
    bench_cmd = ["nice", "-n", "10", str(binary), "--mode", "shard-r27", "--left-size", "4", "--start", "0", "--stop", "1024",
                 "--shard-id", bench["id"], "--out-dir", str(run / "results"), "--source-sha", source_sha, "--checkpoint-every", "128"]
    bench_receipt = invoke(bench_cmd, run, run / "logs" / "benchmark.stdout", run / "logs" / "benchmark.stderr", env)
    bench_checkpoint = run / "results" / checkpoint_name(bench)
    bench_receipt["checkpoint"] = str(bench_checkpoint)
    atomic(run / "benchmark.receipt.json", bench_receipt)
    if bench_receipt["exit_code"] or not bench_checkpoint.exists() or json.loads(bench_checkpoint.read_text()).get("status") != "completed":
        if bench_checkpoint.exists() and json.loads(bench_checkpoint.read_text()).get("status") == "counterexample":
            direct_verify(binary, run, source_sha, bench, bench_checkpoint, env)
        raise SystemExit("benchmark failed/counterexample; production not launched")
    pending, active, halt_left, outcomes = task_list(), {}, set(), []
    while pending or active:
        while len(active) < WORKERS and pending:
            task = pending.popleft()
            if task["left"] in halt_left:
                continue
            cmd = ["nice", "-n", "10", str(binary), "--mode", "shard-r27", "--left-size", str(task["left"]),
                   "--start", str(task["start"]), "--stop", str(task["stop"]), "--shard-id", task["id"],
                   "--out-dir", str(run / "results"), "--source-sha", source_sha, "--checkpoint-every", "256"]
            stdout, stderr = run / "logs" / f"{task['id']}.stdout", run / "logs" / f"{task['id']}.stderr"
            started, t0 = utc(), time.perf_counter()
            out, err = stdout.open("wb"), stderr.open("wb")
            proc = subprocess.Popen(cmd, cwd=run, stdout=out, stderr=err, env=env)
            active[proc.pid] = {"proc": proc, "task": task, "out": out, "err": err, "stdout": stdout, "stderr": stderr,
                                "utc_started": started, "t0": t0}
        finished=[]
        for pid, item in active.items():
            rc = item["proc"].poll()
            if rc is not None: finished.append((pid, rc, item))
        if not finished:
            time.sleep(0.10); continue
        for pid, rc, item in finished:
            del active[pid]; item["out"].close(); item["err"].close(); task=item["task"]
            checkpoint=run / "results" / checkpoint_name(task)
            receipt={"pid":pid,"task":task,"command":["nice","-n","10",str(binary)],"exit_code":rc,
                     "utc_started":item["utc_started"],"utc_finished":utc(),"wall_seconds":time.perf_counter()-item["t0"],
                     "stdout_sha256":sha(item["stdout"]),"stderr_sha256":sha(item["stderr"]),"checkpoint":str(checkpoint)}
            if checkpoint.exists(): receipt["checkpoint_sha256"]=sha(checkpoint); receipt["checkpoint_status"]=json.loads(checkpoint.read_text()).get("status")
            atomic(run / "logs" / f"{task['id']}.receipt.json", receipt); outcomes.append(receipt)
            if rc == 3 and checkpoint.exists() and receipt.get("checkpoint_status") == "counterexample":
                halt_left.add(task["left"])
                direct_verify(binary, run, source_sha, task, checkpoint, env)
                # Stop only active workers for this affected family; their last atomic receipt remains evidence.
                for other in active.values():
                    if other["task"]["left"] == task["left"]: other["proc"].terminate()
            elif rc != 0 or receipt.get("checkpoint_status") != "completed":
                halt_left.add(task["left"])
    atomic(run / "production.outcomes.json", {"source_sha256":source_sha,"halted_left_sizes":sorted(halt_left),"outcomes":outcomes})
    if halt_left:
        raise SystemExit(f"production halted for left sizes {sorted(halt_left)}; do not aggregate")
    aggregate_cmd=[sys.executable, str(ROOT / "resume_aggregate.py"), "aggregate", "--results", str(run / "results"),
                   "--source-sha", source_sha, "--out", str(run / "strict-aggregate.json")]
    aggregate=invoke(aggregate_cmd, run, run / "aggregate.stdout", run / "aggregate.stderr", env)
    atomic(run / "aggregate.receipt.json", aggregate)
    if aggregate["exit_code"]: raise SystemExit("strict coverage aggregation failed")
    atomic(run / "completed.json", {"utc_finished":utc(),"source_sha256":source_sha,"binary_sha256":sha(binary),"status":"completed-strict-coverage"})


if __name__ == "__main__":
    main()
