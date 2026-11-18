"""Rebuild a published area and all common-pose selections from the frozen inputs."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


def replay(config, reference, full):
    comp = Path(__file__).resolve().parent
    expected = json.loads(reference.read_bytes())
    assert sha256(config.read_bytes()).hexdigest() == expected["config_sha256"]
    paths = []
    for name, checksum in expected["placement"]["certificate_sha256"].items():
        assert Path(name).name == name
        path = reference.parent / name
        assert sha256(path.read_bytes()).hexdigest() == checksum
        paths.append(path)
    command = [sys.executable, "-IB", comp/"certify_sigma154_deeper_composition.py",
               comp, config, *paths]
    if full is not None:
        command.extend(["--replay-inputs", full])
    actual = json.loads(subprocess.check_output(command))
    assert actual["placement"]["replayed"] == (full is not None)
    # The extra independent run only changes this provenance flag.
    expected["placement"]["replayed"] = full is not None
    assert actual == expected
    for case in actual["geometry"]["cases"]:
        print(case["label"], case["universal_cover"], case["area"]["decimal"][1])


if __name__ == "__main__":
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--replay-inputs", type=Path)
    args = parser.parse_args()
    replay(args.config, args.reference, args.replay_inputs)
