#!/usr/bin/env python3
"""Integer verification of exterior connectivity of a two-bit bent tensor."""
import argparse
from collections import Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path
import platform
import random
import struct
import time

import check_components as base

SEED = 810828
SAMPLE_COUNTS = {4:64, 5:4}


def datasets(n,exhaustive):
    if exhaustive:
        return range(1 << (1 << n))
    rng = random.Random(SEED+n)
    return [rng.getrandbits(1<<n) for _ in range(SAMPLE_COUNTS[n])]


def truth_product(n,truth):
    # The added sign is (-1)^(z1*z2); the last block is complemented.
    size = 1 << n
    return sum((truth ^ ((1 << size)-1 if z==3 else 0)) << (z*size)
               for z in range(4))


def tensor_support(n,truth):
    size = 1 << n
    support = base.walsh_mask(n,truth)
    expected = sum(support << (size*z) for z in range(4))
    actual = base.walsh_mask(n+2,truth_product(n,truth))
    if expected != actual:
        raise AssertionError('support differs from true Walsh computation')
    return support,actual


def all_kernel_connectivity(n,support,parent):
    new_n = n+2
    count = Counter()
    signature = support & 0xffffffff
    first = None
    for codim,H in base.spaces(new_n):
        exterior = parent & ~H
        components = base.small_components(new_n,exterior) if new_n <= 5 else base.compute_components(new_n,exterior)
        k = len(components)
        if support & (support-1) and k > 1:
            raise AssertionError((n,parent,H,components,'noncharacter bent tensor disconnected'))
        if k>1:
            dim = new_n-codim
            if first is None or dim < first['dimH']:
                first = dict(dimH=dim,kernel_mask=H,exterior_components=components)
        count[new_n-codim,k] += 1
        signature = (signature*1664525 + k + 31*codim) & 0xffffffff
    return count,first,signature


def verify_main(n,exhaustive):
    truths = datasets(n,exhaustive)
    counts = Counter()
    firsts = {}
    cache = {}
    checksum = hashlib.sha256()
    singleton_tables = 0
    for truth in truths:
        support,parent = tensor_support(n,truth)
        if support not in cache:
            cache[support] = all_kernel_connectivity(n,support,parent)
        hist,first,digest = cache[support]
        counts.update(hist)
        if first is not None:
            firsts.setdefault(first['dimH'],dict(base_truth=truth,base_support_mask=support,**first))
        singleton_tables += not bool(support & (support-1))
        checksum.update(struct.pack('<I',digest))
    return dict(base_dimension=n,ambient_dimension=n+2,exhaustive=exhaustive,
                truth_tables=len(truths),singleton_base_tables=singleton_tables,
                kernels=len(base.spaces(n+2)),
                table_kernel_pairs=len(truths)*len(base.spaces(n+2)),
                noncharacter_table_kernel_pairs=(len(truths)-singleton_tables)*len(base.spaces(n+2)),
                distinct_base_supports=len(cache),
                first_disconnect_by_dim=firsts,
                histogram=[dict(dimH=h,components=k,count=c) for (h,k),c in sorted(counts.items())],
                signature_sha256=checksum.hexdigest())


def verify_deletions():
    n=2
    base_support,parent = tensor_support(n,1)
    if base_support != 15:
        raise AssertionError('chosen base is not a two-variable bent phase')
    histogram=Counter()
    witness={}
    for codim,H in base.spaces(n+2):
        dimH=n+2-codim
        exterior=parent&~H
        for d in (0,1,2,3):
            for deleted in combinations(base.points(exterior),d):
                remaining=exterior^sum(1<<u for u in deleted)
                comps=base.small_components(n+2,remaining)
                histogram[dimH,d,len(comps)] += 1
                if len(comps)>1:
                    if d<=2:
                        raise AssertionError((H,deleted,comps,'fewer than three deletions disconnect'))
                    witness.setdefault((dimH,d),dict(kernel_mask=H,removed=deleted,
                                                     components=[list(base.points(c)) for c in comps]))
    return dict(base_truth_mask=1,ambient_dimension=n+2,
                table_kernel_deletion_cases=sum(histogram.values()),
                counts=[dict(dimH=h,n_removed=d,components=k,count=c)
                        for (h,d,k),c in sorted(histogram.items())],
                disconnect_witnesses=[dict(dimH=k[0],n_removed=k[1],**v) for k,v in sorted(witness.items())])


def bundle(components_path):
    print('import types,sys')
    print('module=types.ModuleType("check_components")')
    print('sys.modules["check_components"]=module')
    print('exec(compile(' + repr(Path(components_path).read_text()) + ',"check_components","exec"),module.__dict__)')
    print(Path(__file__).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('exhaustive','sample','deletions','bundle'))
    parser.add_argument('components_path',nargs='?')
    args=parser.parse_args()
    if args.mode=='bundle':
        bundle(args.components_path)
        return
    wall,cpu = time.perf_counter(),time.process_time()
    if args.mode=='deletions':
        records=verify_deletions()
    else:
        dims = (1,2,3) if args.mode=='exhaustive' else SAMPLE_COUNTS
        records=[verify_main(n,args.mode=='exhaustive') for n in dims]
    wall,cpu = time.perf_counter()-wall,time.process_time()-cpu
    print(json.dumps(dict(status='pass',host=platform.node(),python=platform.python_version(),mode=args.mode,
                          records=records,wall_seconds=wall,process_cpu_seconds=cpu,
                          mean_cpu_percent_of_one_core=100*cpu/wall),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
