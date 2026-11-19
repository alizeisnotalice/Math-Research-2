"""Check the two concrete five-body placements, the fixed symmetry host and the perimeter estimate."""
from argparse import ArgumentParser
from pathlib import Path
from hashlib import sha256
import json
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
PACKAGES = {
    "fivefold.zip": "a7586bf42ee989c889184c0e50480d37b4a0fb8649edf268e40c94b4d0ff92c3",
    "basins.zip": "94f11e6c6cf13fd8b56761c43f208569b396133a60f922fc605431cf2279a23a",
    "perimeter.zip": "a89d788b523849e1e964176bb9d784d21c0f4ad8074d63fe55c91d8ac14f7805",
    "symmetric-host.zip": "1fd7030027b2b3943e52c716fdf527d657a37f39c92f1510696d666ac45b9d1b",
    "five-obstacle.zip": "7658ccd50ffbd17b76d029fa62acdd3a1564aea7aee0324f044144f2736f2802",
    "two-poses.zip": "292f63bc7bda4a34e9bb9d539b4a3c2fef9073351ba56dc8f6523b18854057ad",
    "seven-symmetry.zip": "bbd33f54390dd27c9d075581af24a088e519fa69ad12718705c9c3384d89a37a"
}


def run(folder, script, *params, isolated=True):
    flags = ["-I", "-B"] if isolated else ["-B"]
    return subprocess.run([sys.executable, *flags, script, *map(str, params)], cwd=folder,
                          check=True, capture_output=True, text=True).stdout


def recover(folder, script, receipt):
    output = run(folder, script)
    assert output.encode() == (folder / receipt).read_bytes(), script
    return json.loads(output)


def main():
    ap = ArgumentParser(description=__doc__)
    for name in ["tree", "cache", "work"]:
        ap.add_argument(name, type=Path)
    args = ap.parse_args()
    assert __debug__
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    for name, expected in PACKAGES.items():
        archive = HERE / name
        assert sha256(archive.read_bytes()).hexdigest() == expected
        with ZipFile(archive) as z:
            assert all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in z.namelist())
            z.extractall(work / Path(name).stem)
    run(work / "basins/r286-common-basin", "verify.py", HERE.parent,
        args.tree.resolve(), work / "basins-replay", args.cache.resolve())
    basin = json.loads((work / "basins-replay/receipt.json").read_bytes())
    assert basin["basin"]["p30_closed_neighborhood"]["integer_verified_cells"] == 8192

    perimeter = work / "perimeter"
    frozen = {name: (perimeter / name).read_bytes() for name in ["perimeter-result.json", "containment-result.json"]}
    run(perimeter, "containment.py", isolated=False)  # Imports and recomputes the perimeter first.
    assert all(data == (perimeter / name).read_bytes() for name, data in frozen.items())
    geometry = json.loads(frozen["perimeter-result.json"])
    obstruction = recover(work / "five-obstacle", "area_verify.py", "five-receipt.json")
    common = recover(work / "symmetric-host", "verify.py", "receipt.json")
    seven = recover(work / "seven-symmetry", "verify.py", "receipt.json")
    choices = recover(work / "two-poses", "multiple.py", "receipt.json")
    fivefold = recover(work / "fivefold", "verify.py", "receipt.json")
    assert fivefold["boxes"] == 302 and fivefold["balanced_weights"] == 22
    assert all(entry["area"] == common["area"] for entry in [seven, choices, fivefold])
    assert seven["rows"] == 8120 and seven["mixed_tests"] == 1
    assert choices["classes"][0]["paired_tests"] == 395
    # The extension packages can store different auxiliary fields alongside the identical 60 support numbers.
    caps = []
    for name in ["symmetric-host", "seven-symmetry", "two-poses", "fivefold"]:
        data = json.loads((work / name / "host.json").read_bytes())
        caps.append(data["cap"])
    assert caps[0] == caps[1] == caps[2] and len(caps[0]) == 60
    result = {
        "two_fixed_positions": [basin["fixed_six"]["leaf"], basin["basin"]["leaf"]],
        "fixed_six_W_lower": ".838247", "fixed_six_W_strict_upper": ".839503799377581",
        "second_area_upper_with_W_and_five_points": ".837324828897371",
        "second_P30_angular_cells": 8192, "symmetry_host_area": common["area"],
        "cyclic_seven_joint_margin": seven["minimum"], "two_pose_classes": choices["classes"],
        "five_point_obstruction": obstruction, "fivefold_continuum_obstruction": fivefold,
        "perimeter_upper": geometry["perimeter_bound"], "density_upper": geometry["q"],
        "blend_upper_area": geometry["area_at_t"],
    }
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
