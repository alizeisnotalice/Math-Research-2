"""Compare Gaussian quadrature orders at identical binary32 inputs."""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify import f32, integrate


def load(path):
    return json.loads(Path(path).read_text())


def sign(x):
    return int(x > 0) - int(x < 0)


def main():
    source = load("points.json")
    reports = load("matched-gpu.json")
    values = {r["order"]: r["values"] for r in reports}
    cases = []
    margin = 2e-6
    scaled_keys = ("exact","leading")
    for i, item in enumerate(source):
        p = {k: (v if k=="tuned" else f32(v)) for k,v in item["point"].items()}
        cpu = integrate(p)
        direct = {k:cpu[k] for k in scaled_keys}
        gpu = {order:{k:f32(calc[i][k]) for k in scaled_keys} for order,calc in values.items()}
        cases.append({
            **item, "cpu":direct, "gpu":gpu,
            "separated_cpu": min(abs(x) for x in direct.values()) >= margin,
            "gpu_error": {
                order:max(abs(val[k]-direct[k]) for k in scaled_keys)
                for order,val in gpu.items()
            },
            "both_signs_agree": {
                order: all(sign(val[k])==sign(direct[k]) for k in scaled_keys)
                for order,val in gpu.items()
            },
        })
    summary={
        "count":len(cases),
        "tested_orders":list(values),
        "native_command_seconds":{r["order"]:r["commandSeconds"] for r in reports},
        "max_error_vs_cpu":{
            order:max(row["gpu_error"][order] for row in cases)
            for order in values
        },
        "max_order_difference":{
            f"{left}-{right}":max(max(abs(row["gpu"][left][k]-row["gpu"][right][k]) for k in scaled_keys)
                                   for row in cases)
            for left,right in ((96,128),(128,160))
        },
        "separated_cpu":sum(row["separated_cpu"] for row in cases),
        "separated_cpu_sign_inconsistencies":{
            order:sum(row["separated_cpu"] and not row["both_signs_agree"][order] for row in cases)
            for order in values
        },
        "cases":cases,
    }
    Path("matched-summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print({k:v for k,v in summary.items() if k!="cases"})


if __name__=="__main__":
    main()
