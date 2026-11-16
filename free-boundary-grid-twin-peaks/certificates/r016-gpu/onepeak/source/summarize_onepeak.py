#!/usr/bin/env python3
"""Summarize saved one-peak SMC repeats; SE is descriptive only."""
import glob, json, math, os, statistics, sys

rows=[]
for path in sorted(glob.glob(sys.argv[1])):
    data=[json.loads(x) for x in open(path)]
    reps=[x for x in data if x["kind"]=="replicate"]
    z=[x["Z"] for x in reps]
    rows.append({"file":os.path.basename(path),"n":data[0]["n"],"P":data[0]["particles"],"R":len(reps),
                 "seed_domain_hex":data[0]["seed_domain_hex"],"mean_Z":statistics.mean(z),
                 "sample_variance_Z":statistics.variance(z) if len(z)>1 else 0.0,
                 "descriptive_SE_Z":statistics.stdev(z)/math.sqrt(len(z)) if len(z)>1 else 0.0,
                 "dead":sum(x["dead"] for x in reps),"wall_seconds":data[-1]["wall_seconds"],
                 "gpu_seconds_sum":sum(x["gpu_hardware_seconds"] for x in reps),
                 "synchronous_dispatch_seconds_sum":sum(x["synchronous_dispatch_seconds"] for x in reps)})
for x in rows: print(json.dumps(x,sort_keys=True))
