#!/usr/bin/env python3
"""Bounded parallel exact search followed by independent full-ledger checking.

Input is a list of {tag, min_total, total, max_terms}. Each first degree is a
separate process. No result status is inferred from timeout or process failure.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('config', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--seconds', type=float, default=14400)
    args = ap.parse_args()
    assert Path('/Volumes/data').is_mount()
    assert 1 <= args.workers and 0 < args.seconds <= 14400
    configs = json.loads(args.config.read_text())
    assert configs and len({c['tag'] for c in configs}) == len(configs)
    for c in configs:
        assert c['tag'].isalnum() and 2 <= c['min_total'] <= c['total']
        assert c['max_terms'] >= 2
    args.out.mkdir(parents=True, exist_ok=True)
    base = Path(__file__).resolve().parent
    start = time.monotonic()
    deadline = start + args.seconds

    def run(task):
        c, m = task
        out = args.out / f"{c['tag']}-m{m:02d}"
        out.mkdir()  # Refuse accidental overwriting or duplicate live runs.
        state = {**c, 'm': m, 'status': 'started', 'process_exited': False}
        commands = [
            [sys.executable, str(base/'support_search_layer.py'),
             '--total', str(c['total']), '--min-total', str(c['min_total']),
             '--min-degree', str(m), '--max-first-degree', str(m),
             '--max-terms', str(c['max_terms']),
             '--seconds', str(max(1, deadline-time.monotonic())), '--out', str(out)],
            [sys.executable, str(base/'verify_support_layer.py'), str(out)]]
        for phase, command in zip(('search', 'verify'), commands):
            if time.monotonic() >= deadline:
                state.update(status='unknown_batch_deadline', phase=phase)
                break
            with (out/(phase+'.stdout')).open('w') as stdout, (out/(phase+'.stderr')).open('w') as stderr:
                child = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                    env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1',
                         'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
                state.update(phase=phase, pid=child.pid, status='running', process_exited=False)
                (out/'job.json').write_text(json.dumps(state))
                try:
                    code = child.wait(timeout=max(0.01, deadline-time.monotonic()))
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
                    state.update(status='unknown_timeout', process_exited=True)
                    break
                state['process_exited'] = True
                if code:
                    state.update(status='unknown_failure', returncode=code)
                    break
        else:
            state['status'] = 'verified'
        (out/'job.json').write_text(json.dumps(state, indent=2)+'\n')
        print(json.dumps(state), flush=True)
        return state

    tasks = [(c,m) for c in configs for m in range(1, c['total']//2+1)]
    def cost(task):
        c,m = task
        def count(d):
            return sum(math.comb(d-1,j) for j in range(min(d-1,c['max_terms']-2)+1))
        return sum(count(m)*count(t-m) for t in range(c['min_total'],c['total']+1) if t >= 2*m)
    tasks.sort(key=cost, reverse=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run,task) for task in tasks]
        jobs = [f.result() for f in as_completed(futures)]
    report = {'workers':args.workers, 'configs':configs, 'jobs':jobs,
              'elapsed_seconds':time.monotonic()-start,
              'all_workers_exited':True}
    (args.out/'batch.json').write_text(json.dumps(report,indent=2)+'\n')
    if any(j['status'] != 'verified' for j in jobs):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
