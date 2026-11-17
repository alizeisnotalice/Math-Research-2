#!/usr/bin/env python3
"""Bounded remote driver for the two authorized R14 D18 capacity cases.

The driver runs foreground child processes, records their exit status, and
samples only its own DP executable plus host CPU idleness.  It deliberately
does not enumerate or launch any candidate other than the two named below.
"""
import datetime as dt
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
SOURCE = ROOT / "capacity_dp_parallel.cpp"
INPUT = ROOT / "eight-point-candidates.json"
BIN = ROOT / "capacity_dp_parallel"
OUT = ROOT / "d18-run-summary.json"
NAMES = ("remote_d18_free_pair", "remote_d18_mixed")
THREADS = 15
LIMIT_SECONDS = 900
SAMPLE_SECONDS = 60
VALIDATIONS = (
    {
        "name": "old_serial_triangular_Q4_D7",
        "D": 7,
        "expected_kappa": 8,
        "basis": [
            "0xf000f000f000f", "0xf000f000f000f0", "0xf000f000f000f00",
            "0xf000f000f000f000", "0xaaaaaaaaaaaaaaaa", "0xcccccccccccccccc",
            "0xc0c0c0c0c0c0c0c0",
        ],
        "reference": "round013/capacity-verification.json triangular_Q4",
    },
    {
        "name": "old_serial_selector_task11_D13",
        "D": 13,
        "expected_kappa": 16,
        "basis": [
            12297829382473034410, 7392283487217560234, 3512861310540840960,
            868082073854803968, 287332900645175040, 67779170971733196,
            4279075918001100, 224561496588288, 56913480056832, 4042322160,
            817901760, 202116108, 50529027,
        ],
        "reference": "round014/selector-capacity-results.jsonl id=11",
    },
)


def sha256(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def command_output(argv):
    p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       text=True, check=False)
    return {"argv": argv, "returncode": p.returncode,
            "stdout": p.stdout, "stderr": p.stderr}


def cpu_sample(label, child_pid=None):
    top = command_output(["/usr/bin/top", "-l", "1", "-n", "0"])
    cpu_line = next((x.strip() for x in top["stdout"].splitlines()
                     if "CPU usage:" in x), "")
    ps = command_output(["/bin/ps", "-axo", "pid,ppid,pcpu,pmem,etime,command"])
    rows = [x for x in ps["stdout"].splitlines()
            if str(BIN) in x and "grep" not in x]
    own_cpu = 0.0
    for row in rows:
        parts = row.split(None, 5)
        if len(parts) >= 3:
            try:
                own_cpu += float(parts[2])
            except ValueError:
                pass
    return {"utc": dt.datetime.now(dt.timezone.utc).isoformat(), "label": label,
            "child_pid": child_pid, "host_cpu_line": cpu_line,
            "own_dp_cpu_percent_sum": own_cpu, "own_dp_ps_rows": rows}


def run_case(spec, mode):
    command = [str(BIN), str(THREADS), str(spec["D"])] + [str(x) for x in spec["basis"]]
    began = time.monotonic()
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    rec = {"name": spec["name"], "mode": mode, "D": spec["D"],
           "expected_kappa": spec.get("expected_kappa"), "reference": spec.get("reference"),
           "candidate_upper": spec.get("candidate_upper"), "command": command,
           "started_utc": started, "timeout_seconds": LIMIT_SECONDS, "cpu_samples": []}
    proc = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    rec["pid"] = proc.pid
    rec["cpu_samples"].append(cpu_sample("start", proc.pid))
    timed_out = False
    next_sample = time.monotonic() + SAMPLE_SECONDS
    while proc.poll() is None:
        now = time.monotonic()
        if now - began >= LIMIT_SECONDS:
            timed_out = True
            proc.terminate()
            try:
                stdout, stderr = proc.communicate(timeout=15)
            except subprocess.TimeoutExpired:
                proc.kill()
                stdout, stderr = proc.communicate()
            break
        if now >= next_sample:
            rec["cpu_samples"].append(cpu_sample("interval", proc.pid))
            next_sample += SAMPLE_SECONDS
        time.sleep(min(2.0, max(0.05, next_sample - now)))
    else:
        stdout, stderr = proc.communicate()
    rec["cpu_samples"].append(cpu_sample("end", proc.pid))
    rec.update({"finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "elapsed_seconds_wall": time.monotonic() - began,
                "returncode": proc.returncode, "timed_out": timed_out,
                "stdout": stdout, "stderr": stderr})
    try:
        rec["parsed"] = json.loads(stdout.strip())
    except json.JSONDecodeError:
        rec["parsed"] = None
    rec["matches_expected_serial_kappa"] = (
        rec["parsed"] is not None and spec.get("expected_kappa") is not None and
        rec["parsed"].get("kappa") == spec["expected_kappa"] and proc.returncode == 0 and not timed_out)
    return rec


def main():
    input_obj = json.loads(INPUT.read_text())
    candidates = input_obj if isinstance(input_obj, list) else input_obj["candidates"]
    chosen = [x for x in candidates if x.get("name") in NAMES]
    if tuple(x.get("name") for x in chosen) != NAMES:
        raise RuntimeError("authorized D18 candidates missing or input order changed")
    if any(x.get("D") != 18 for x in chosen):
        raise RuntimeError("D18 guard failed")
    result = {"schema": "a28-r014-d18-remote-v1", "authorized_names": list(NAMES),
              "threads": THREADS, "per_case_timeout_seconds": LIMIT_SECONDS,
              "source_sha256": sha256(SOURCE), "input_sha256": sha256(INPUT),
              "driver_sha256": sha256(__file__),
              "host": command_output(["/bin/hostname"]),
              "compiler": command_output(["/usr/bin/clang++", "--version"]),
              "python": command_output([sys.executable, "--version"]),
              "started_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "runs": []}
    compile_rec = command_output(["/usr/bin/clang++", "-O3", "-std=c++17", "-pthread",
                                  str(SOURCE), "-o", str(BIN)])
    result["compile"] = compile_rec
    if compile_rec["returncode"] != 0:
        result["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return 20
    for spec in VALIDATIONS:
        result["runs"].append(run_case(spec, "validation"))
        OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not all(x["matches_expected_serial_kappa"] for x in result["runs"]):
        result["validation_passed"] = False
        result["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return 21
    result["validation_passed"] = True
    for candidate in chosen:
        result["runs"].append(run_case(candidate, "authorized_d18"))
        OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    result["finished_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    result["all_authorized_runs_ok"] = all(
        x["returncode"] == 0 and not x["timed_out"] and x["parsed"] is not None
        for x in result["runs"] if x["mode"] == "authorized_d18")
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0 if result["all_authorized_runs_ok"] else 22


if __name__ == "__main__":
    sys.exit(main())
