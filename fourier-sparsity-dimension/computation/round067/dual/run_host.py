"""Run one all-thread CPU shard and concurrent independently evaluated Metal shard."""
import argparse
import csv
import hashlib
import json
import pathlib
import re
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sample_gpu():
    try:
        q=subprocess.run(["/usr/sbin/ioreg","-r","-c","IOAccelerator","-d","1"],capture_output=True,text=True,timeout=3)
    except subprocess.TimeoutExpired:
        return {"epoch":time.time(),"ioreg_exit_code":None,"metrics":{"Device Utilization %":[]}}
    fields={}
    for name in ("Device Utilization %","Renderer Utilization %","Tiler Utilization %"):
        fields[name]=[int(x) for x in re.findall(re.escape('"'+name+'"')+r"\s*=\s*([0-9]+)",q.stdout)]
    return {"epoch":time.time(),"ioreg_exit_code":q.returncode,"metrics":fields}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("out",type=pathlib.Path)
    p.add_argument("first",type=int)
    p.add_argument("count",type=int)
    args=p.parse_args()
    out=args.out.resolve()
    out.mkdir(exist_ok=True)
    root=pathlib.Path(__file__).resolve().parent
    cores=int(subprocess.check_output(["/usr/sbin/sysctl","-n","hw.logicalcpu"],text=True))
    c,metal=[root/name for name in ("clique_cpu","clique_metal")]
    hashes={name:sha(root/name) for name in ("clique_cpu.c","AuditCommon.h","ClickKernel.metal","CliqueMetal.swift","run_host.py")}
    gate=out/"start_epoch.txt"
    if gate.exists(): raise RuntimeError("Old gate")
    (out/"ready.json").write_text(json.dumps(dict(hashes=hashes,logical_cores=cores,time=time.time()))+"\n")
    deadline=time.monotonic()+180
    while not gate.exists():
        if time.monotonic()>deadline: raise RuntimeError("Gate not issued")
        time.sleep(.05)
    epoch=float(gate.read_text())
    before=sample_gpu()
    while time.time()<epoch: time.sleep(min(.01,max(0,epoch-time.time())))
    begun=time.time()
    cpu_json=out/"cpu.json"
    cpu_csv=out/"cpu_usage.csv"
    gpu_json=out/"gpu.json"
    arguments=[str(args.first),str(args.count)]
    with (out/"work.log").open("w") as log:
        cpu=subprocess.Popen([str(c)]+arguments+[str(cores),str(epoch),str(cpu_json),str(cpu_csv)],stdout=log,stderr=subprocess.STDOUT)
        gpu=subprocess.Popen([str(metal),str(gpu_json)]+arguments+["2048"],cwd=root,stdout=log,stderr=subprocess.STDOUT)
        print("R67 PIDs",cpu.pid,gpu.pid,"CPU threads",cores,flush=True)
        gpu_monitor=[]
        while cpu.poll() is None or gpu.poll() is None:
            gpu_monitor.append(dict(sample_gpu(),cpu_running=cpu.poll() is None,gpu_running=gpu.poll() is None))
            time.sleep(.14)
        result_codes=[cpu.wait(),gpu.wait()]
    ended=time.time()
    after=sample_gpu()
    existing=subprocess.run(["/bin/ps","-p",f"{cpu.pid},{gpu.pid}","-o","pid="],capture_output=True,text=True).stdout.split()
    result={
        "schema":"a28-r067-dual-run-receipt-v1","same_requested_epoch":epoch,"launched_epoch":begun,
        "ended_epoch":ended,"logical_cores":cores,"cpu_pid":cpu.pid,"metal_pid":gpu.pid,
        "cpu_gpu_exit_codes":result_codes,"remaining_task_process_pids":existing,"source_sha256":hashes,
        "cpu_binary_sha256":sha(c),"metal_binary_sha256":sha(metal),
        "input_interval_half_open":[args.first,args.first+args.count],
        "GPU_telemetry_scope":"Whole device, includes unrelated work",
        "GPU_before":before,"GPU_after":after,"GPU_samples":gpu_monitor,
    }
    if result_codes==[0,0]:
        a=json.loads(cpu_json.read_bytes())
        b=json.loads(gpu_json.read_bytes())
        shared = ("first_id","count","strict_graph","strict_basis","strict_activity",
                  "sums","positional_checksum64")
        assert all(a[key] == b[key] for key in shared)
        metrics=list(csv.DictReader(cpu_csv.read_text().splitlines()))
        workers=[s for s in metrics if int(s["active_workers"])==cores]
        if not workers:workers=metrics
        result.update(GPU_CPU_aggregate_exact_match=True,
                      cpu_gpu_outputs_sha256=[sha(cpu_json),sha(gpu_json)],
                      sampled_cpu_telemetry_sha256=sha(cpu_csv),
                      peak_a28_cpu_process_percent=max(float(r["a28_process_cpu_percent"]) for r in workers),
                      peak_whole_machine_cpu_percent=max(float(r["whole_machine_cpu_active_percent"]) for r in workers),
                      peak_whole_device_gpu_percent=max((n for sample in gpu_monitor if sample["gpu_running"]
                                                        for n in sample["metrics"]["Device Utilization %"]),default=None))
    with (out/"receipt.json").open("w") as file: json.dump(result,file,sort_keys=True,indent=2)
    print({key:value for key,value in result.items() if key!="GPU_samples"},flush=True)
    if result_codes!=[0,0] or existing: raise SystemExit(1)


if __name__=="__main__":
    main()
