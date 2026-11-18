"""Check one representative per signing modulo addition of a linear form."""
import json
import multiprocessing as mp
import platform
import socket
import sys
import time
from signed_flags import check

def batch(bounds):
    first, last, side = bounds
    started = time.process_time()
    supports = reps = signs = ties = 0
    for i in range(first, last):
        assigned_to_mini = i % 3 == 0
        if assigned_to_mini != (side == 'mini'):
            continue
        mask = i << 1
        span = {0}
        free = mask
        for e in range(1, 16):
            if mask >> e & 1 and e not in span:
                free ^= 1 << e
                span |= {j ^ e for j in span}
        supports += 1
        sub = free
        while True:
            count, _ = check(mask, sub)
            reps += 1
            signs += len(span)
            ties += count
            if not sub:
                break
            sub = (sub - 1) & free
    return supports, reps, signs, ties, time.process_time() - started

if __name__ == '__main__':
    workers, side = int(sys.argv[1]), sys.argv[2]
    assert side in ('book', 'mini')
    start = time.perf_counter()
    jobs = [(s, min(s+16,32768), side) for s in range(1,32768,16)]
    with mp.get_context('fork').Pool(workers) as pool:
        totals = list(map(sum, zip(*pool.imap_unordered(batch,jobs,chunksize=1))))
    print(json.dumps(dict(host=socket.gethostname(), python=platform.python_version(),
        workers=workers, side=side, supports=totals[0], representatives=totals[1],
        signings=totals[2], tied_cuts=totals[3], cpu_seconds=totals[4],
        wall_seconds=time.perf_counter()-start), indent=2))
