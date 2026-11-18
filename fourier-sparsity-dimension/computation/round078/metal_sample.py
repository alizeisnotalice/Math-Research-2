#!/usr/bin/env python3
"""Cross-check exact selector spectra with direct integer Metal transforms."""
import argparse
from array import array
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import time
import types


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError('cannot load ' + name)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def stream_code(audit_path, metal_path):
    for name, path in [('audit', audit_path), ('metal_backend', metal_path)]:
        source = Path(path).read_text()
        print('import types, sys')
        print(f'module = types.ModuleType({name!r})')
        print(f'sys.modules[{name!r}] = module')
        print(f'exec(compile({source!r}, {name!r}, "exec"), module.__dict__)')
    print(Path(__file__).read_text().rsplit('\nif __name__ ==',1)[0])
    print('run_check()', end='\n')


def run_check():
    from audit import indicator_bands, r49_cells, weak_folding_cells
    from metal_backend import Metal
    records = []
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    metal = Metal()
    try:
        for name, constructor in [('R49',r49_cells),('Tarannikov',weak_folding_cells)]:
            for height in (1,2):
                core, cells = constructor(height)
                n = core + len(cells)
                location = {x:i for i,cell in enumerate(cells) for x in cell}
                n_values = 1 << n
                coefficients, scale, _ = indicator_bands(cells,core)
                # Pack physical truth values independently from the spectral indicator expansion.
                truth = array('I')
                for word in range(0,n_values,32):
                    packed_word = 0
                    for bit in range(32):
                        point = word + bit
                        tile = location[point % (1 << core)]
                        packed_word |= ((point >> (core + tile)) & 1) << bit
                    truth.append(packed_word)
                raw, gpu_seconds, wait_seconds = metal.compute(truth,n_values,1)
                actual = array('i')
                actual.frombytes(raw)
                if len(actual) != n_values:
                    raise AssertionError('output length differs')
                for index, w in enumerate(actual):
                    if w * scale != coefficients.get(index,0) * n_values // 2:
                        raise AssertionError((name,height,index,w,coefficients.get(index,0)))
                records.append(dict(name=f'{name}-{height}',dimension=n,walsh_coefficients_checked=n_values,
                                    gpu_command_seconds=gpu_seconds,commit_and_wait_seconds=wait_seconds,
                                    int32_walsh_half_sha256=hashlib.sha256(raw).hexdigest()))
    finally:
        metal.close()
    wall_seconds = time.perf_counter()-wall_start
    gpu = sum(r['gpu_command_seconds'] for r in records)
    wait = sum(r['commit_and_wait_seconds'] for r in records)
    print(json.dumps(dict(status='pass',python=platform.python_version(),host=platform.node(),
                          device=metal.device_name,backend='Metal int32 Walsh sums',records=records,
                          wall_seconds=wall_seconds,process_cpu_seconds=time.process_time()-cpu_start,
                          total_gpu_command_seconds=gpu,gpu_command_fraction_of_wall=gpu/wall_seconds,
                          gpu_command_fraction_of_commit_and_wait=gpu/wait),indent=2,sort_keys=True))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle',action='store_true',help='emit a self-contained script to stdout')
    parser.add_argument('audit_path')
    parser.add_argument('metal_backend_path')
    args=parser.parse_args()
    if args.bundle:
        stream_code(args.audit_path,args.metal_backend_path)
    else:
        load_module('audit',args.audit_path)
        load_module('metal_backend',args.metal_backend_path)
        run_check()
