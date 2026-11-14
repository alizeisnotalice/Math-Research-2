#!/usr/bin/env python3
"""Supply portable paths to the byte-preserved R002 search programs."""
import argparse
import importlib.util
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('kind', choices=('directional', 'position'))
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline', type=Path,
                        default=Path(__file__).parent / 'inputs/legacy_r004_transfer_34_27.json')
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    sys.path.insert(0, str(directory))
    name = 'scan_perturbations' if args.kind == 'directional' else 'scan_position_grid'
    spec = importlib.util.spec_from_file_location(name, directory / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    args.output.mkdir(parents=True, exist_ok=True)
    module.OUT = args.output.resolve()
    if args.kind == 'directional':
        module.INPUT = args.input.resolve()
    else:
        module.BEST = args.input.resolve()
        module.BASE = args.baseline.resolve()
    return module.main()


if __name__ == '__main__':
    sys.exit(main())
