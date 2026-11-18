"""Check lexicographic signed cuts against every rank-deficit inequality in F2^4."""
import argparse
import hashlib
import json
import multiprocessing as mp
import platform
import random
import socket
import time

def subspaces():
    spaces = {1}
    for v in range(1, 16):
        for h in list(spaces):
            spaces.add(h | sum(1 << (i ^ v) for i in range(16) if h >> i & 1))
    return sorted(spaces)

SPACES = subspaces()
DIM = {h: h.bit_count().bit_length() - 1 for h in SPACES}
CUTS = [sum(1 << e for e in range(1, 16) if (e & v).bit_count() % 2)
        for v in range(16)]
EXT = {h: [k for k in SPACES if k != h and k & h == h] for h in SPACES}

def check(mask, signs):
    h = 1
    candidates = list(range(16))
    layers = []
    while mask & ~h:
        best = None
        for k in EXT[h]:
            size = (mask & k & ~h).bit_count()
            rank = DIM[k] - DIM[h]
            if best is None or size * best[1] > best[0] * rank or (
                size * best[1] == best[0] * rank and rank > best[1]
            ):
                best = size, rank, k
        size, rank, k = best
        layer = mask & k & ~h
        values = {v: ((CUTS[v] ^ signs) & layer).bit_count() for v in candidates}
        maximum = max(values.values())
        candidates = [v for v in candidates if values[v] == maximum]
        layers.append((h, k, layer, size, rank))
        h = k
    for v in candidates:
        x = (CUTS[v] ^ signs) & mask
        for h, k, layer, size, rank in layers:
            for t in [h] + EXT[h]:
                if t & k == t:
                    count = (x & layer & ~t).bit_count()
                    assert 2 * count >= (layer & ~t).bit_count()
                    assert 2 * rank * count >= size * (DIM[k] - DIM[t])
    return len(candidates), len(layers)

def batch(args):
    start, end, salt = args
    begin = time.process_time()
    ties = depth = 0
    digest = hashlib.sha256()
    for i in range(start, end):
        signs = random.Random((salt << 16) + i).getrandbits(16)
        m, n = check(i << 1, signs)
        ties += m
        depth = max(depth, n)
        digest.update(bytes((m, n)))
    return ties, depth, digest.hexdigest(), time.process_time() - begin

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('workers', type=int)
    p.add_argument('salt', type=int)
    args = p.parse_args()
    jobs = [(s, min(s+128,32768),args.salt) for s in range(1,32768,128)]
    start = time.perf_counter()
    with mp.get_context('fork').Pool(args.workers) as pool:
        results = pool.map(batch, jobs, chunksize=1)
    print(json.dumps(dict(host=socket.gethostname(), python=platform.python_version(),
        workers=args.workers, salt=args.salt, supports=32767,
        max_depth=max(x[1] for x in results), tied_cuts=sum(x[0] for x in results),
        digest=hashlib.sha256(''.join(x[2] for x in results).encode()).hexdigest(),
        cpu_seconds=sum(x[3] for x in results), wall_seconds=time.perf_counter()-start),indent=2))
