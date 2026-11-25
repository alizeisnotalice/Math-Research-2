"""Replay the archived rounds without writing files. Run with Python -B."""

import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def checked(script, *args):
    return subprocess.check_output(
        [sys.executable, "-B", str(HERE / script), *args], cwd=HERE
    )


if __name__ == "__main__":
    subprocess.run(
        [
            sys.executable,
            "-B",
            str(HERE.parent / "replay.py"),
            "--check",
            str(HERE.parent / "certificate.json"),
        ],
        check=True,
    )
    for name, script, args in (
        ("kernel", "lattice", ["3", "8", "32"]),
        ("short", "short", []),
        ("minors", "minors", []),
        ("signs", "signs", []),
        ("platform", "platform", []),
        ("apt", "valuations", ["apt"]),
        ("mixed", "valuations", ["mixed"]),
        ("riv", "valuations", ["riv"]),
        ("cmf", "cmf", []),
        ("gap", "gap", []),
    ):
        if json.loads(checked(script + ".py", *args)) != json.loads(
            (HERE / (name + ".json")).read_text()
        ):
            raise ValueError(name + " certificate mismatch")
        print(name, "OK", flush=True)
    checked("check32.py")
    print("Independent kernel verifier OK")
