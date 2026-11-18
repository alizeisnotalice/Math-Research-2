"""Compile and run one exact phase-corrected flat orbit shard, with a common UTC start gate."""
import argparse
import hashlib
import json
import os
import pathlib
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("first", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("swift_cache", type=pathlib.Path)
    parser.add_argument("clang_cache", type=pathlib.Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(exist_ok=True)
    temp = out / "temp"
    temp.mkdir(exist_ok=True)
    gate = out / "start_epoch.txt"
    if gate.exists():
        raise RuntimeError("Start gate already exists")
    exe = out / "r53_flat_orbit_metal"
    source = pathlib.Path(__file__).resolve().with_name("R53FlatOrbitMetal.swift")
    env = os.environ.copy()
    env.update({key: str(temp) for key in ("TMPDIR","TMP","TEMP","XDG_CACHE_HOME")})
    env["CLANG_MODULE_CACHE_PATH"] = str(args.clang_cache)
    compile_args = ["/usr/bin/swiftc","-O","-module-cache-path",str(args.swift_cache),str(source),"-o",str(exe)]
    compile_started = time.monotonic()
    compilation = subprocess.run(compile_args,capture_output=True, text=True,env=env)
    compile_seconds = time.monotonic()-compile_started
    print(compilation.stderr,end="",flush=True)
    if compilation.returncode != 0:
        raise SystemExit(compilation.returncode)
    ready = {"ready_unix":time.time(),"source_sha256":digest(source),"binary_sha256":digest(exe)}
    (out/"ready.json").write_text(json.dumps(ready)+"\n")
    stage_end=time.monotonic()+80
    while not gate.is_file():
        if time.monotonic()>stage_end:
            raise RuntimeError("Start gate not released")
        time.sleep(0.03)
    target = float(gate.read_text())
    while time.time()<target:
        time.sleep(min(.015,max(0,target-time.time())))
    wall_start=time.time()
    file=out/"result.json"
    seed=source.with_name("r11_minimal_input.json")
    command=[str(exe),str(seed),str(file),str(args.first),str(args.end)]
    with (out/"stdout.txt").open("w") as log:
        job=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env)
        print("R53 worker",job.pid,"range",args.first,args.end,flush=True)
        exit_code=job.wait()
    wall_end=time.time()
    active=subprocess.run(["/bin/ps","-p",str(job.pid),"-o","pid="],capture_output=True,text=True)
    version=subprocess.run(["/usr/bin/swiftc","--version"],capture_output=True,text=True)
    result={
        "schema":"a28-r053-flat-orbit-dual-shard-receipt-v1","start_epoch":target,
        "actual_start_unix":wall_start,"actual_end_unix":wall_end,
        "compile_exit_code":compilation.returncode,"compile_seconds":compile_seconds,
        "swift_version":(version.stderr+version.stdout).strip(),
        "run_exit_code":exit_code,"pid":job.pid,"exited":not active.stdout.strip(),
        "source_sha256":ready["source_sha256"],"binary_sha256":ready["binary_sha256"],
        "runner_sha256":digest(pathlib.Path(__file__)),
        "data_sha256":digest(file) if file.exists() else None,
        "run_arguments":[str(args.first),str(args.end)],"input_sha256":digest(seed),
        "worker_stdout":(out/"stdout.txt").read_text(),
    }
    (out/"receipt.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,sort_keys=True))
    if exit_code != 0 or not result["exited"]:
        raise SystemExit(exit_code or 1)


if __name__ == "__main__":
    main()
