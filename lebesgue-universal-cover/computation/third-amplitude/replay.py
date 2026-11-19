"""Recompute the three Fourier-amplitude hosts and the two analytic fixed hosts."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import zipfile

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
DIGESTS = {
    "low": "19c9638c1a7f837cab3036a6f9147b2c62d215d5d38b16ffb3adaaaa1fc4c1b5",
    "middle": "1f32c8271bd726851539772c2d1e0fc8d0200ad3a1485d173b9c78f7c6d1f79f",
    "high": "5136b3598230eb8d1f138a8a555d0812064d93626fd1de0ae09d2e9b18284808",
    "stability": "6022b7d789e2ebd8d10357c62bc8f3e373f51e5e5d145a716a1ed75f38637adc",
    "same-host": "0552a8aa3b9c229c1ae2f1c66fa59d404df401e90cb784a5c553b31056e816d7",
}


def expand(name, work):
    package = HERE / f"{name}.zip"
    assert hashlib.sha256(package.read_bytes()).hexdigest() == DIGESTS[name]
    destination = work / name
    with zipfile.ZipFile(package) as bundle:
        for member in bundle.namelist():
            path = Path(member)
            assert not path.is_absolute() and ".." not in path.parts
        bundle.extractall(destination)
    return destination


def run(folder, script, *args):
    return subprocess.run(
        [sys.executable, "-I", "-B", script, *map(str, args)],
        cwd=folder, check=True, capture_output=True, text=True
    ).stdout


def repeat(folder, script, record):
    output = run(folder, script)
    assert output.encode() == (folder / record).read_bytes(), record
    return json.loads(output)


def analytic_hosts():
    pi = arb.pi()
    root3 = arb(3).sqrt()
    budget = arb("0.844082046521218581802948589770767704354")
    coefficient = arb(3023) / 10000
    t = arb(4264) / 100000
    delta = (2 / root3 - 1) * t
    area = pi / 8 + root3 / 4 + (pi / 4 - root3 / 2) * t + (pi / 8 - root3 / 4) * t**2
    perimeter = pi / 2 + root3 + (pi / 2 - root3) * t
    tube = area + perimeter * delta + pi * delta**2
    assert coefficient * arb(".02182") < delta and tube < budget
    # The isoperimetric support inequality gives P(U) <= 4 area(U).
    upper = arb("0.844088795387818568068160578195374139951")
    t_blend = arb(9892) / 10000
    blend = upper - (1 - t_blend)**2 * (upper - pi / 4)
    q_limit = (1 - ((upper - budget) / (upper - pi / 4)).sqrt()) / 2
    assert t_blend / 2 < q_limit and blend < budget
    return {
        "tube_radius": str(delta), "tube_area": str(tube), "tube_margin": str(budget - tube),
        "tube_q_floor_branch": ".02132", "tube_deficit": ".02182",
        "universal_blend_q": ".4946", "universal_blend_area_upper_enclosure": str(blend),
        "universal_blend_q_budget_crossing": str(q_limit)
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="new output directory")
    args = parser.parse_args()
    assert __debug__, "use a Python interpreter with assertions enabled"
    args.work.mkdir(parents=True, exist_ok=False)
    work = args.work.resolve()
    ctx.prec = 192
    folders = {name: expand(name, work) for name in DIGESTS}
    low_dir = folders["low"] / "green-farkas-low-face-hosts"
    run(low_dir, "verify_release.py", work / "low.json")
    low = json.loads((work / "low.json").read_bytes())["continuous_hosts"]["directed60r7435.json"]
    middle = repeat(folders["middle"], "verify.py", "receipt.json")
    miss = repeat(folders["middle"], "inradius.py", "inradius-receipt.json")
    high = repeat(folders["high"], "final_verify.py", "final-receipt.json")
    linear = repeat(folders["stability"], "adaptive.py", "receipt.json")
    joint_dir = folders["same-host"] / "green-r3-merge"
    run(joint_dir, "verify.py", work / "joint.json")
    joint = json.loads((work / "joint.json").read_bytes())
    assert joint["status"] == "PASS" and joint["stability"] == linear
    assert joint["old"]["r7435"]["certified_to_exceed_B_if_retaining_old_host"]
    assert arb(miss["difference"]) > arb(".00000650194")
    assert Fraction(low["third_radius"]) == Fraction(".07435")
    assert arb(middle["lower_critical"]) < arb(".073945")
    assert arb(middle["upper_critical"]) > arb(".163066")
    assert arb(high["dual_critical"]) < arb(".16259")
    assert Fraction(".073945") < Fraction(".07435") < Fraction(".16259") < Fraction(".163066")
    assert all(row.get("triangles", row.get("complete_triples")) == 8120 for row in (low, middle, high))
    assert all(row.get("pairs", row.get("opposite_pairs")) == 30 for row in (low, middle, high))
    areas = [arb(row["area"]) for row in (low, middle, high)]
    upper_three = arb("0.844076541772819091385619585133383")
    assert all(area < upper_three for area in areas)
    assert upper_three < arb("0.844082046521218581802948589770767704354")
    result = {"three_host_areas": list(map(str, areas)), "three_host_strict_upper": str(upper_three),
              "angle_ranges": [["0", ".07435"], [".073945", ".163066"], [".16259", "2/pi"]],
              "same_low_high_host_is_a_separate_construct": True, **analytic_hosts()}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
