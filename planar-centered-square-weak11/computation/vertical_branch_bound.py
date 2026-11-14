#!/usr/bin/env python3
"""Single-process certified branch and bound for the vertical scaling family.

The complete leaf partition, including pruned leaves, is written as JSON.
All decisions and bounds are rational; elapsed time is used only to stop.
"""
import argparse
import hashlib
import heapq
import json
import os
import platform
import sys
import time
from fractions import Fraction as F
from pathlib import Path

from vertical_scaling import VerticalFamily


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temp = path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,indent=2)+'\n')
    os.replace(temp,path)


def prepare(data_path, cache_path):
    input_hash = digest(data_path)
    family_hash = digest(Path(__file__).with_name('vertical_scaling.py'))
    fraction_fields = ('bmin','bmax','radius','px','py','lam','mass','tail_bound')
    if cache_path and Path(cache_path).exists():
        data = json.loads(Path(cache_path).read_text())
        if data['input_sha256'] != input_hash or data['family_code_sha256'] != family_hash:
            raise ValueError('cache input or family code hash mismatch')
        f = VerticalFamily.__new__(VerticalFamily)
        for name in fraction_fields:
            setattr(f,name,F(data['attributes'][name]))
        for name in ('scale','bits','enumeration'):
            setattr(f,name,data['attributes'][name])
        f.boxes = [tuple(box) for box in data['boxes']]
        return f, 'cache'
    f = VerticalFamily(json.loads(Path(data_path).read_text()))
    if cache_path:
        attrs = {name:str(getattr(f,name)) for name in fraction_fields}
        attrs.update({name:getattr(f,name) for name in ('scale','bits','enumeration')})
        atomic_json(cache_path, {'input_sha256':input_hash,'family_code_sha256':family_hash,'attributes':attrs,'boxes':f.boxes})
    return f, 'enumeration'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--output', required=True)
    parser.add_argument('--target-gap', default='1/100000000')
    parser.add_argument('--seconds', type=float, default=30)
    parser.add_argument('--checkpoint-seconds', type=float, default=30)
    parser.add_argument('--cache', help='JSON cache created by this script')
    parser.add_argument('--resume', action='store_true', help='resume from --output checkpoint')
    args = parser.parse_args()
    gap = F(args.target_gap)
    if gap <= 0 or args.seconds <= 0 or args.checkpoint_seconds <= 0:
        parser.error('positive target gap and time limit required')
    start = time.monotonic()
    family, preparation = prepare(args.input,args.cache)
    init_seconds = time.monotonic()-start
    factor = family.lam/family.mass
    samples = {}
    incumbent, best_beta = None, None

    def sample(beta):
        nonlocal incumbent, best_beta
        if beta not in samples:
            lower, upper = family.evaluate(beta)
            samples[beta] = (lower*factor, upper*factor)
            if incumbent is None or lower*factor > incumbent:
                incumbent, best_beta = lower*factor, beta

    low, high = family.bmin, family.bmax
    heap, finished = [], []
    sequence = divisions = 0
    if args.resume:
        previous = json.loads(Path(args.output).read_text())
        if previous['input_sha256'] != digest(args.input) or previous['code_sha256'] != digest(__file__) or previous['family_code_sha256'] != digest(Path(__file__).with_name('vertical_scaling.py')):
            raise ValueError('checkpoint input or code hash mismatch')
        if F(previous['target_gap']) != gap:
            raise ValueError('checkpoint target gap mismatch')
        samples = {F(b):(F(v['lower']),F(v['upper'])) for b,v in previous['samples'].items()}
        incumbent, best_beta = F(previous['incumbent_lower']), F(previous['best_beta'])
        divisions = previous['divisions']
        for leaf in previous['leaves']:
            a,b,u = map(F,(leaf['lo'],leaf['hi'],leaf['upper']))
            if leaf['status'] == 'pruned':
                finished.append((a,b,u))
            elif leaf['status'] == 'active':
                sequence += 1
                heapq.heappush(heap,(-u,sequence,a,b))
            else:
                raise ValueError('invalid checkpoint leaf status')
    else:
        sample((low+high)/2)
        upper = family.interval_upper(low, high)*factor
        heap = [(-upper, 0, low, high)]

    def checkpoint(reason):
        leaves = [(a,b,u,'pruned') for a,b,u in finished]
        leaves += [(a,b,-neg_u,'active') for neg_u,_,a,b in heap]
        leaves.sort()
        assert leaves[0][0] == low and leaves[-1][1] == high
        assert all(leaves[i][1] == leaves[i+1][0] for i in range(len(leaves)-1))
        global_upper = max([incumbent]+[u for _,_,u,_ in leaves])
        assert all(u < incumbent for _,_,u,status in leaves if status == 'pruned')
        if reason == 'target_gap' and global_upper-incumbent > gap:
            raise AssertionError('target gap is not established')
        out = {'input_sha256':digest(args.input), 'code_sha256':digest(__file__), 'family_code_sha256':digest(Path(__file__).with_name('vertical_scaling.py')), 'python':sys.version, 'platform':platform.platform(), 'beta_min':str(low), 'beta_max':str(high), 'target_gap':str(gap), 'seconds_limit':args.seconds, 'checkpoint_seconds':args.checkpoint_seconds, 'preparation':preparation, 'initialization_seconds':init_seconds, 'elapsed_seconds':time.monotonic()-start, 'bits':family.bits, 'tail_bound':str(family.tail_bound), 'enumeration':family.enumeration, 'status':reason, 'divisions':divisions, 'best_beta':str(best_beta), 'incumbent_lower':str(incumbent), 'incumbent_lower_float':float(incumbent), 'global_upper':str(global_upper), 'global_upper_float':float(global_upper), 'certified_gap':str(global_upper-incumbent), 'certified_gap_float':float(global_upper-incumbent), 'samples':{str(b):{'lower':str(v[0]),'upper':str(v[1])} for b,v in sorted(samples.items())}, 'leaves':[{'lo':str(a),'hi':str(b),'upper':str(u),'status':status} for a,b,u,status in leaves]}
        atomic_json(args.output,out)
        print(json.dumps({k:v for k,v in out.items() if k not in ('samples','leaves')},indent=2),flush=True)

    reason = 'time_limit'
    next_checkpoint = start+args.checkpoint_seconds
    while heap:
        # Pruning remains valid even if the incumbent later improves.
        while heap and -heap[0][0] < incumbent:
            neg_u, _, a, b = heapq.heappop(heap)
            finished.append((a,b,-neg_u))
        if not heap:
            reason = 'target_gap'
            break
        if -heap[0][0]-incumbent <= gap:
            reason = 'target_gap'
            break
        if time.monotonic()-start >= args.seconds:
            break
        neg_u, _, a, b = heapq.heappop(heap)
        mid = (a+b)/2
        # Sample the parent midpoint; the first midpoint is already cached.
        sample(mid)
        for left,right in ((a,mid),(mid,b)):
            bound = family.interval_upper(left,right)*factor
            sequence += 1
            if bound < incumbent:
                finished.append((left,right,bound))
            else:
                heapq.heappush(heap,(-bound,sequence,left,right))
        divisions += 1
        if time.monotonic() >= next_checkpoint:
            checkpoint('running')
            next_checkpoint = time.monotonic()+args.checkpoint_seconds
    checkpoint(reason)


if __name__ == '__main__':
    main()
