"""R75 CPU execution, paired Cantor resonances and four-device timing."""
from datetime import datetime,timedelta
import json, math, re, statistics, sys
from pathlib import Path

import mpmath as mp


def dt(v):
    return datetime.fromisoformat(v)


def main():
    mp.mp.dps=42
    cpu=json.loads(Path(sys.argv[1]).read_text())
    mac_cpu=json.loads(Path(sys.argv[2]).read_text())
    clock=json.loads(Path(sys.argv[3]).read_text())
    local_native=Path(sys.argv[4]).read_text()
    mac_life=Path(sys.argv[5]).read_text()
    mac_json=json.loads(Path(sys.argv[6]).read_text())
    mac_ps=Path(sys.argv[7]).read_text()
    coarse=json.loads(Path(sys.argv[8]).read_text())
    cases=[c for shard in cpu["shards"] for c in shard["tests"]]
    rows=[r for c in cases for r in c["rows"]]
    worker_pids={s["pid"] for s in cpu["shards"]}
    lo,hi=clock["post_run_remote_minus_local_ns"]
    local_start=max(dt(s["start"]) for s in cpu["shards"])
    local_end=min(dt(s["finish"]) for s in cpu["shards"])
    mac_cpu_start=dt(mac_cpu["all_worker_timers_start_upper"])-timedelta(microseconds=lo//1000)
    mac_cpu_end=dt(mac_cpu["all_worker_timers_end_lower"])-timedelta(microseconds=-(-hi//1000))
    local_gpu_start,local_gpu_pid=re.search(r"(?m)^start=(\S+).*pid=(\d+)",local_native).groups()
    local_gpu_end,local_gpu_total=re.search(r"(?m)^end=(\S+).*gpu_seconds=(\S+)",local_native).groups()
    assert re.search(r"(?m)^exit=0 pid="+local_gpu_pid+"$",local_native)
    mac_gpu_start,mac_gpu_pid,mac_gpu_end,exit_code=re.search(
        r"(?m)^start (\S+)\npid (\d+)\nend (\S+) exit (\d+)$",mac_life).groups()
    assert exit_code=="0"
    positives=[]
    for s in re.split(r"(?m)^sample ",mac_ps)[1:]:
        if re.search(r"(?m)^\s*"+mac_gpu_pid+r"\s",s):
            positives.append(dt(s.splitlines()[0]))
    # The first process sample is between its marker and next marker, which
    # in this particular evidence still leaves ample lead time to the window.
    mac_bound_start=positives[1]+timedelta(seconds=1,microseconds=-(lo//1000))
    mac_bound_end=positives[-1]-timedelta(microseconds=-(-hi//1000))
    first=max(local_start,mac_cpu_start,dt(local_gpu_start),mac_bound_start)
    last=min(local_end,mac_cpu_end,dt(local_gpu_end),mac_bound_end)
    all_obs=[];common=[]
    for item in cpu["samples"]:
        workers=[x for x in item["ps"] if int(x[0]) in worker_pids]
        if len(workers)==len(worker_pids) and all(x[4][0]=="R" for x in workers):
            p=(dt(item["utc"]),sum(float(x[2]) for x in workers),item["machine_cpu_pct"][2])
            all_obs.append(p)
            if first<=p[0]<=last:
                common.append(p)
    candidates=[(dt(t),float(s)) for t,s in re.findall(r"(?m)^utc=(\S+).*gpu_seconds=(\S+)",local_native)]
    candidates.append((dt(local_gpu_end),float(local_gpu_total)))
    inside=[x for t,x in candidates if first<=t<=last]
    local_command=inside[-1]-inside[0]
    mac_wall_upper=(dt(mac_gpu_end)+timedelta(seconds=1)-dt(mac_gpu_start)).total_seconds()
    mac_commands=float(mac_json["gpuSeconds"])
    mac_slack=mac_wall_upper-mac_commands

    by_alpha={}
    for alpha in sorted({c["alpha"] for c in cases}):
        cc=[c for c in cases if c["alpha"]==alpha]
        biggest=[(c,*[next(r for r in reversed(c["rows"]) if r["base_harmonic"]==k)
                        for k in (1,next(j["base_harmonic"] for j in c["rows"] if j["base_harmonic"]>1))])
                 for c in cc]
        ratios=[]
        for c,suppressed,_ in biggest:
            A=c["digits"];n=suppressed["n"];m=suppressed["exponent"]
            gap_count=sum(e>d+1 for d,e in zip(A,A[1:]))
            h=(2*mp.pi**2*mp.mpf(c["c"])*n*n)**(-1/mp.mpf(alpha))
            predicted=2*gap_count*h*mp.gamma(1+1/mp.mpf(alpha))*(len(A)**m-1)/(len(A)-1)
            ratios.append(float(mp.mpf(suppressed["coefficient"][0])/predicted))
        by_alpha[alpha]={"models":len(cc),"suppressed_deepest_ratio_to_coarse_asymptotic_median":
                         statistics.median(ratios),
                         "surviving_deepest_phase_ratio_error_median":statistics.median(
                                 float(r["phase_ratio_error"]) for _,_,r in biggest)}
    answer={
      "hardware":"Apple M6, 12 logical cores","platform":cpu["platform"],"python":cpu["python"],
      "mpmath":cpu["mpmath"],"scipy":cpu["scipy"],"cpu_pid":cpu["parent"],
      "cpu_worker_pids":sorted(worker_pids),"cpu_utc":[cpu["start"],cpu["finish"]],
      "cpu_worker_timer_utc":[local_start.isoformat(),local_end.isoformat()],
      "cpu_wall_s":cpu["wall_s"],"worker_cpu_s":sum(r["cpu_seconds"] for r in cpu["shards"]),
      "cpu_worker_exit_verified":cpu["worker_exit_verified"],
      "all_cpu_worker_runnable_samples":len(all_obs),"cpu_runnable_pct_median":statistics.median(x[1] for x in all_obs),
      "system_cpu_idle_range_in_all_worker_samples":[min(x[2] for x in all_obs),max(x[2] for x in all_obs)],
      "cases":len(cases),"frequencies":len(rows),"passed":sum(r["cross_pass"] for r in rows),
      "max_cross_abs":max(r["cross_difference"] for r in rows),
      "max_omitted_levels_bound":max(r["tail_bound"] for r in rows),
      "independent_coarse_cases":len(coarse["rows"]),
      "coarse_max_abs_cross":max(float(r["abs_difference"]) for r in coarse["rows"]),
      "coarse_max_rel_cross":max(float(r["relative_difference"]) for r in coarse["rows"]),
      "clock_remote_minus_local_bounds_ns_from_mac_cpu_receipt":[lo,hi],
      "clock_assumption":"Post-run bounds treated as stable through the workloads.",
      "four_job_process_overlap_utc":[first.isoformat(),last.isoformat()],
      "four_job_process_overlap_s":(last-first).total_seconds(),
      "all_12_local_cpu_samples_in_shared_window":len(common),
      "local_gpu_commands_inside_shared_window_s_lower":local_command,
      "mac_gpu_sequential_command_time_s":mac_commands,
      "mac_gpu_lifecycle_s_upper":mac_wall_upper,
      "dual_gpu_command_time_lower_from_global_totals_s":max(0,local_command-mac_slack),
      "mac_gpu_pid_positive_marker_utc_range":[positives[0].isoformat(),positives[-1].isoformat()],
      "mac_gpu_pid_positive_markers":len(positives),
      "by_alpha":by_alpha}
    if len(sys.argv)==11:
        mac_cycles=Path(sys.argv[9]).read_text()
        mtime_ns=int(Path(sys.argv[10]).read_text().strip())
        zero=dt("1970-01-01T00:00:00+00:00")
        mtime_upper=zero+timedelta(microseconds=-(-mtime_ns//1000))
        init_lower=dt(mac_gpu_start)-timedelta(microseconds=-(-hi//1000))
        init_upper=(mtime_upper-timedelta(seconds=float(mac_json["wallSeconds"]),
                        microseconds=lo//1000)+timedelta(microseconds=2))
        cycles=[(0,0.,0.)]+[
            (int(j),float(g),float(w)) for j,g,w in re.findall(
            r"(?m)^cycle=(\d+).*gpu_seconds=(\S+) wall=(\S+)$",mac_cycles)]
        settled=[]
        for (j0,g0,w0),(j1,g1,w1) in zip(cycles,cycles[1:]):
            cycle_lower=init_lower+timedelta(seconds=w0,microseconds=-2)
            cycle_upper=init_upper+timedelta(seconds=w1,microseconds=2)
            if first<=cycle_lower<=cycle_upper<=last:
                settled.append((j1,g1-g0,cycle_lower.isoformat(),cycle_upper.isoformat()))
        mac_inside=math.fsum(c[1] for c in settled)
        answer.update({
            "mac_cycle_timer_origin_local_bounds":[init_lower.isoformat(),init_upper.isoformat()],
            "mac_cycle_indices_guaranteed_inside":[settled[0][0],settled[-1][0]],
            "mac_cycle_commands_inside_shared_window_s_lower":mac_inside,
            "dual_gpu_command_time_inside_shared_window_s_lower":
              max(0.,local_command+mac_inside-(last-first).total_seconds()),
        })
    Path("r075-local-cpu-receipt.json").write_text(json.dumps(answer,indent=2)+"\n")
    print(json.dumps({k:v for k,v in answer.items() if k!="by_alpha"},indent=2))


if __name__=="__main__":
    main()
