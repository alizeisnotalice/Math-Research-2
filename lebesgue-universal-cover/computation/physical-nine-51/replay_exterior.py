"""Replay the .83751 exterior from the 1,024 seed parts and the native layer receipts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile


HERE = Path(__file__).resolve().parent


def unpack(name, destination):
    with ZipFile(HERE / name) as archive:
        for member in archive.namelist():
            path = Path(member)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError(member)
        archive.extractall(destination)


def run(script, *args):
    subprocess.run([sys.executable, "-B", str(script), *map(str, args)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("computation", type=Path)
    parser.add_argument("tree", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    comp = args.computation.resolve()
    old_core = comp / "global-083750" / "records" / "retained.json"
    program = out / "outside"
    raw = out / "outer238"
    unpack("outside-512.zip", program)
    unpack("outside-238.zip", raw)
    seed = out / "seed"
    run(program / "seal_seeds.py", "unpack", seed, HERE / "seeds.zip")
    layer_sets = [out / group for group in ("layers-a", "layers-b", "layers-c")]
    for folder in layer_sets:
        unpack(folder.name + ".zip", folder)
    with (out / "root.json").open("w") as log:
        subprocess.run([sys.executable, "-B", str(program / "audit_root_high.py"),
                        str(comp), str(args.tree.resolve())], check=True, stdout=log)
    assert (out / "root.json").read_bytes() == (program / "root-083751-083755.json").read_bytes()
    layers = out / "layers"
    run(program / "audit51_replays.py", seed, old_core, layers, *layer_sets)
    for name in ("retained.json", "audit.json"):
        assert (layers / name).read_bytes() == (program / "stage4" / name).read_bytes()
    info = json.loads((layers / "audit.json").read_bytes())
    assert info["parts"] == 1024 and info["still_in_seed"] == 5954
    group = program / "groups01"
    direct = [out / "group01.json"]
    run(program / "audit_new_outer.py", group / "index.json", direct[0],
        group / "five-083751.jsonl", group / "five-b-083751.jsonl", group / "qonly-b-083751.jsonl")
    for name in ("group2", "group3"):
        folder = raw / name
        dest = out / (name + ".json")
        run(raw / "audit_new_outer.py", folder / "index.json", dest,
            folder / "five-083751.jsonl", folder / "qonly-083751.jsonl")
        direct.append(dest)
    result = out / "exterior.json"
    run(program / "complete_51_exterior.py", old_core, layers, result, *direct)
    assert result.read_bytes() == (program / "exterior-083751-complete.json").read_bytes()


if __name__ == "__main__":
    main()
