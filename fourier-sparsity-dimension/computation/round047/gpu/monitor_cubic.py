"""Observe system GPU use while a prepared cubic Walsh Metal sweep runs."""
import argparse
import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess
import time


def read_gpu():
    results = subprocess.run(["/usr/sbin/ioreg", "-r", "-c", "IOAccelerator", "-d", "1"],
                             text=True, capture_output=True, timeout=2)
    metrics = {}
    for name in ("Device Utilization %", "Renderer Utilization %", "Tiler Utilization %"):
        metrics[name] = [int(i) for i in re.findall(re.escape('"' + name + '"') + r"\s*=\s*([0-9]+)",
                                                   results.stdout)]
    return {"epoch_unix": time.time(), "measurements": metrics, "ioreg_exit_code": results.returncode}


def sha(file):
    return hashlib.sha256(pathlib.Path(file).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("first", type=int)
    parser.add_argument("count", type=int)
    parser.add_argument("seed")
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--sample-stride", type=int, default=32)
    opts = parser.parse_args()
    directory = opts.directory.absolute()
    directory.mkdir(exist_ok=True)
    gate = directory / "start_epoch.txt"
    if gate.exists():
        raise RuntimeError("This sweep needs a new output directory")
    root = pathlib.Path(__file__).resolve().parent
    binary = root / "r47_cubic_walsh_gpu"
    source = root / "R47CubicWalshLoad.swift"
    assert (directory / "pilot.json").is_file()
    pilot = json.loads((directory / "pilot.json").read_bytes())
    if not pilot["parseval_and_sample_checks_passed"]:
        raise RuntimeError("Pilot unverified")
    ready = {"source_sha256": sha(source), "binary_sha256": sha(binary), "pilot_sha256": sha(directory / "pilot.json")}
    (directory / "ready.json").write_text(json.dumps(ready) + "\n")
    print("Ready, supervisor pid", os.getpid(), flush=True)
    deadline = time.monotonic() + 75
    while not gate.is_file():
        if time.monotonic() > deadline:
            raise RuntimeError("Synchronized start gate not supplied")
        time.sleep(0.035)
    epoch = float(gate.read_text())
    before = read_gpu()
    while time.time() < epoch:
        time.sleep(min(0.015, max(0, epoch - time.time())))
    result = directory / "cubic_sweep_result.json"
    log = directory / "cubic_sweep.log"
    command = [str(binary), str(result), str(opts.first), str(opts.count), opts.seed,
               str(opts.batch_size), str(opts.sample_stride)]
    start = time.time()
    with log.open("w") as out:
        job = subprocess.Popen(command, stdout=out, stderr=subprocess.STDOUT)
        print("GPU process", job.pid, "started at UTC Unix", start, flush=True)
        telemetry = []
        while job.poll() is None:
            telemetry.append(read_gpu())
            time.sleep(0.16)
    end = time.time()
    after = read_gpu()
    status = subprocess.run(["/bin/ps","-p",str(job.pid),"-o","pid="],
                            text=True,capture_output=True)
    samples = [n for sample in telemetry for n in sample["measurements"]["Device Utilization %"]]
    receipt = {
        "schema": "a28-r047-cubic-metal-load-receipt-v1",
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scheduled_unix_epoch": epoch,
        "observed_start_epoch": start, "observed_end_epoch": end,
        "working_wall_seconds": end-start,
        "process_id": job.pid, "process_exited": not status.stdout.strip(),
        "exit_code": job.returncode, "source": ready,
        "result_sha256": sha(result) if result.exists() else None,
        "runner_sha256": sha(__file__),
        "reproduction_arguments": [str(opts.first), str(opts.count), opts.seed,
                                   str(opts.batch_size), str(opts.sample_stride)],
        "baseline_ioaccelerator": before, "postrun_ioaccelerator": after,
        "sampled_while_task_alive": telemetry,
        "system_device_utilization_percent_bounds": [min(samples),max(samples)] if samples else None,
        "system_samples_at_least_90_percent": sum(n >= 90 for n in samples),
        "utilization_scope": "Whole device; includes processes other than this scanner.",
    }
    file = directory / "cubic_sweep_receipt.json"
    file.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print({k:v for k,v in receipt.items() if k != "sampled_while_task_alive"}, flush=True)
    if job.returncode or not receipt["process_exited"]:
        raise SystemExit(job.returncode or 1)


if __name__ == "__main__":
    main()
