"""Reduce finite-contact jump results and process telemetry."""
from collections import defaultdict
from decimal import Decimal
import json
from pathlib import Path
import statistics
import sys


def dec(s):
    return Decimal(s)


def run(path):
    data=json.loads(path.read_text())
    cases=[c for s in data["shards"] for c in s["cases"]]
    pids={s["pid"] for s in data["shards"]}
    loads,idle=[],[]
    for r in data["samples"]:
        pool=[p for p in r["ps"] if int(p[0]) in pids]
        if len(pool)==len(pids) and all(p[4].startswith("R") for p in pool):
            loads.append(sum(float(p[2]) for p in pool))
            idle.append(r["whole_cpu_user_sys_idle"][2])
    stats=defaultdict(list)
    for c in cases:
        stats[tuple(c["orders"])].append(c)
    by_order={}
    for orders, group in sorted(stats.items()):
        last=[c["rows"][-1] for c in group]
        by_order[str(orders)]={
            "profiles":len(group),"finest_offset":last[0]["z"],
            "max_finest_jump_error_divided_by_max_1_and_abs_exact_jump": float(max(
                abs(dec(r["jump_error"]))/max(1,abs(dec(c["predicted_jump_divided_by_source"])))
                for r,c in zip(last,group))),
            "max_finest_side_slope_error_divided_by_max_1_and_abs_exact_slope":float(max(
                abs(dec(s["scaled_slope_error"]))/max(1,abs(dec(s["first_slope"])))
                for r in last for s in r["sides"])),
            "sign_matched_at_finest":sum(dec(r["observed_jump"])*dec(c["predicted_jump_divided_by_source"])>0
                       for r,c in zip(last,group)),
        }
    sides=[s for c in cases for r in c["rows"] for s in r["sides"]]
    doc={"start_utc":data["start"],"finish_utc":data["finish"],"parent_pid":data["parent"],
         "wall_seconds":data["wall_s"],"worker_pids":sorted(pids),
         "subprocess_exit_ps_verified":data["worker_exit_verified"],
         "total_worker_cpu_seconds":sum(s["cpu_seconds"] for s in data["shards"]),
         "all_runnable_samples":len(loads),"group_pct_median":statistics.median(loads),
         "group_pct_min_max":[min(loads),max(loads)],"machine_idle_min_max":[min(idle),max(idle)],
         "profile_count":len(cases),"symmetric_pair_quad_tests":len(sides),
         "cross_abs_max":max(s["cross_error"] for s in sides),
         "cross_relative_or_absolute_max":max(s["cross_error"]/max(1,abs(float(s["normalized"]))) for s in sides),
         "by_order":by_order}
    return doc


if __name__=="__main__":
    doc=run(Path(sys.argv[1]))
    text=json.dumps(doc,indent=2)+"\n"
    Path("jump-summary.json").write_text(text)
    print(text,end="")
