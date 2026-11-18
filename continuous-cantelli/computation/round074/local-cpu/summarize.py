"""Compact execution and resonance receipt for the Cantor scan."""
from datetime import datetime
import json, math, statistics, sys
from pathlib import Path


def main():
    d=json.loads(Path(sys.argv[1]).read_text())
    workers={s["pid"] for s in d["shards"]}
    cases=[c for s in d["shards"] for c in s["tests"]]
    rows=[r for c in cases for r in c["rows"]]
    samples=[]
    for s in d["samples"]:
        pool=[p for p in s["ps"] if int(p[0]) in workers]
        if len(pool)==len(workers) and all(p[4].startswith("R") for p in pool):
            samples.append((datetime.fromisoformat(s["utc"]),
                            sum(float(p[2]) for p in pool),s["machine_cpu_pct"][2]))
    by_alpha={}
    for a in sorted({c["alpha"] for c in cases}):
        same=[c for c in cases if c["alpha"]==a]
        by_alpha[a]={
            "number_of_configurations":len(same),
            "median_phase_ratio_error_at_resonant_modes": [
                statistics.median(float(c["rows"][i]["phase_ratio_error"]) for c in same)
                for i in range(3)],
            "median_phase_ratio_error_at_off_resonance":statistics.median(
                float(c["rows"][-1]["phase_ratio_error"]) for c in same),
        }
    out={"cpu_parent_pid":d["parent"],"cpu_worker_pids":sorted(workers),
         "cpu_utc":[d["start"],d["finish"]],"cpu_wall_s":d["wall_s"],
         "worker_exit_verified":d["worker_exit_verified"],
         "worker_cpu_s":sum(s["cpu_seconds"] for s in d["shards"]),
         "all_worker_runnable_samples":len(samples),
         "all_runnable_observed_span_s":(samples[-1][0]-samples[0][0]).total_seconds(),
         "group_median_cpu_pct":statistics.median(s[1] for s in samples),
         "whole_idle_min_max":[min(s[2] for s in samples),max(s[2] for s in samples)],
         "configurations":len(cases),"frequency_rows":len(rows),
         "frequency_rows_passed":sum(r["cross_pass"] for r in rows),
         "maximum_cross_difference":max(r["cross_difference"] for r in rows),
         "maximum_omitted_gap_bound":max(r["tail_bound"] for r in rows),
         "maximum_product_vs_explicit_prefix_error":max(float(r["prefix_four_levels_error"]) for r in rows),
         "minimum_abs_fourier":min(math.hypot(*map(float,r["coefficient"])) for r in rows),
         "by_alpha":by_alpha}
    if len(sys.argv)==3:
        checks=json.loads(Path(sys.argv[2]).read_text())
        out["explicit_physical_gap_crosschecks_max_difference"]=max(
            r["max_absolute_difference"] for r in checks["rows"])
    Path("cantor-summary.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="by_alpha"},indent=2))


if __name__=="__main__":
    main()
