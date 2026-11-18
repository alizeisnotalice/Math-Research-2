"""Cross-check sample values and representative maxima in every Metal stratum."""
import json
from pathlib import Path
from sys import argv

from verify import f32, integrate


def eligible(p, offset):
    tuned = p["tuned"] != 0
    h, b = p["h"], p["b"]
    return {
        "all": True,
        "tuned": tuned,
        "tuned_offset_lt_001": tuned and offset < 0.01,
        "tuned_h_lt_001": tuned and h < 0.01,
        "tuned_h_lt_0001_offset_lt_001": tuned and h < 0.001 and offset < 0.01,
        "untuned": not tuned,
        "tuned_negative_b": tuned and b < 0,
        "tuned_nonnegative_b": tuned and b >= 0,
    }


def sign(x):
    return (x > 0) - (x < 0)


def check(row):
    p = {key: (value if key == "tuned" else f32(value)) for key, value in row["point"].items()}
    raw = {key: f32(row[key]) for key in ("exact", "leading")}
    direct = integrate(p)
    errors = {key: abs(raw[key] - direct[key]) for key in raw}
    away = min(abs(raw["exact"]), abs(raw["leading"])) >= 2e-6
    difference = sign(raw["exact"]) != sign(raw["leading"])
    mode = "near_zero" if not away else ("separated_flip" if difference else "separated_agreement")
    return {
        "index": row["index"], "point": p, "mode": mode, "gpu": raw,
        "cpu": {k: direct[k] for k in raw}, "errors": errors,
        "signs_match": all(sign(raw[key]) == sign(direct[key]) for key in raw),
    }


def main():
    rows = json.loads(Path(argv[1]).read_text())
    buckets = {}
    for row in rows:
        if row["index"] == 984450415:
            offender = check(row)
            continue
        for stratum, include in eligible(row["point"], row["offset"]).items():
            if not include:
                continue
            x, y = row["exact"], row["leading"]
            mode = ("near_zero" if min(abs(x), abs(y)) < 2e-6
                    else "separated_flip" if sign(x) != sign(y) else "separated_agreement")
            modes = buckets.setdefault(stratum, {})
            cases = modes.setdefault(mode, [])
            if len(cases) < 4:
                cases.append(row)

    checks = {}
    for stratum, modes in buckets.items():
        checks[stratum] = {mode: [check(case) for case in cases] for mode, cases in modes.items()}
    total = [case for modes in checks.values() for cases in modes.values() for case in cases]
    separated = [case for case in total if case["mode"] != "near_zero"]
    report = {
        "hardware_rows": len(rows),
        "selected_checks_with_stratum_memberships": len(total),
        "unique_checked_sample_ids": len({case["index"] for case in total} | {offender["index"]}),
        "max_error": max(max(case["errors"].values()) for case in total + [offender]),
        "separated_sign_inconsistencies": sum(not case["signs_match"] for case in separated),
        "strata": checks, "previously_unstable": offender,
    }
    Path(argv[2]).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def rare_main(input_path,output_path):
    rows = json.loads(Path(input_path).read_text())
    checks = [check(row) for row in rows]
    assert checks and all(case["mode"] == "separated_flip" for case in checks)
    out = {
        "checked": len(checks),
        "signs_agree": sum(case["signs_match"] for case in checks),
        "cpu_away_from_margin": sum(min(abs(v) for v in case["cpu"].values()) >= 2e-6
                                   for case in checks),
        "max_error": max(max(case["errors"].values()) for case in checks),
        "minimum_cpu_absolute": min(min(abs(v) for v in case["cpu"].values()) for case in checks),
        "samples": checks,
    }
    Path(output_path).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    if argv[1] == "--rare":
        rare_main(argv[2],argv[3])
    else:
        main()
