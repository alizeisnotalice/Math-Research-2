#!/usr/bin/env python3
"""Exact local producer for the unreceived total-54 cap-six support chunks.

This is deliberately a *producer*, not a proof or an exclusion ledger: it
rebuilds the full plan independently, subtracts only hash-bound persisted
events, and writes the remaining chunks to a new local ledger.
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import gzip
import hashlib
import importlib.util
import itertools
import json
import os
import shutil
import sys
import time
from collections import Counter
from pathlib import Path

OLD = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round043-macbook-support54-corrected/run')
SRC = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/support54_macbook43/plan')
AUDITOR = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/verify_slice45_39/audited_chunk_verifier.py')
OUT_DEFAULT = Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round047-local54-complement')
PLAN_SHA = '29e2d5d8473de4859e4af9d60440519bb2c544a3880fdbfd2b61d090975a1e26'
EVENTS_SHA = '33fa152157a659e63ef44cbf68be7bff4c0352fd82a495255d05835e2ac8facf'
AUDIT_SHA = 'c92536ecc17f8abe9df2dbff26dedb553d8e8e97b18bb48773bf941c06a3f509'
FROZEN_SHA = 'd68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710'
CAP = 6


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def event_sha(event: dict) -> str:
    return hashlib.sha256(json.dumps(event, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def supports(m: int):
    for size in range(2, min(m + 1, CAP) + 1):
        for inside in itertools.combinations(range(1, m), size - 2):
            yield (0, *inside, m)


def independent_full54(plan: dict) -> dict[tuple[int, int, int, tuple[int, ...]], dict]:
    """Generate total54 independently and require literal equality with plan."""
    expected = []
    for m in range(1, 54 // 2 + 1):
        n = 54 - m
        for index, I in enumerate(supports(m)):
            expected.append({'total': 54, 'm': m, 'n': n, 'index': index, 'I': list(I)})
    listed = [j for j in plan['jobs'] if j['total'] == 54]
    if listed != expected:
        raise RuntimeError('total54 plan is not the independently regenerated enumeration')
    ans = {}
    for job in expected:
        key = (job['total'], job['m'], job['n'], job['index'], tuple(job['I']))
        if key in ans:
            raise RuntimeError('independent enumeration duplicate')
        ans[key] = job
    if len(ans) != 101583:
        raise RuntimeError(f'unexpected total54 full count {len(ans)}')
    return ans


def load_received(full: dict[tuple[int, int, int, tuple[int, ...]], dict]) -> set[tuple[int, int, int, tuple[int, ...]]]:
    """Accept only known protocol events and known full-plan chunk keys."""
    allowed_nonchunks = {'support54_local_started', 'support54_remote_started', 'support54_remote_workers'}
    received: set[tuple[int, int, int, tuple[int, ...]]] = set()
    total53 = 0
    for line_no, line in enumerate((OLD / 'events.jsonl').read_text().splitlines(), 1):
        event = json.loads(line)
        kind = event.get('event')
        if kind == 'chunk':
            try:
                key = (int(event['total']), int(event['m']), int(event['n']), int(event['index']), tuple(event['I']))
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError(f'malformed historical chunk event at {line_no}') from exc
            if key[0] == 54:
                if key not in full:
                    raise RuntimeError(f'historical unknown total54 key at {line_no}: {key}')
                if key in received:
                    raise RuntimeError(f'historical duplicate total54 key at {line_no}: {key}')
                received.add(key)
            elif key[0] == 53:
                total53 += 1
            else:
                raise RuntimeError(f'historical chunk outside plan totals at {line_no}')
        elif kind not in allowed_nonchunks:
            raise RuntimeError(f'unknown historical event kind at {line_no}: {kind!r}')
    if total53 != 83681 or len(received) != 10784:
        raise RuntimeError(f'historical counts differ: total53={total53}, total54={len(received)}')
    return received


def atomic_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f'.{os.getpid()}.partial')
    with temp.open('w') as f:
        f.write(json.dumps(obj, sort_keys=True, separators=(',', ':')) + '\n')
        f.flush(); os.fsync(f.fileno())
    os.replace(temp, path)


def append_event(events: Path, event: dict) -> None:
    events.parent.mkdir(parents=True, exist_ok=True)
    with events.open('a') as f:
        f.write(json.dumps(event, sort_keys=True, separators=(',', ':')) + '\n')
        f.flush(); os.fsync(f.fileno())


_CLASSIFIER = None


def classifier(frozen: str):
    global _CLASSIFIER
    if _CLASSIFIER is None:
        spec = importlib.util.spec_from_file_location('frozen54', frozen)
        module = importlib.util.module_from_spec(spec)
        assert spec and spec.loader
        spec.loader.exec_module(module)
        _CLASSIFIER = module
    return _CLASSIFIER


def produce_one(root_s: str, frozen_s: str, manifest_sha: str, job: dict) -> dict:
    """Worker: atomically write a gzip plus its receipt, then return its event."""
    root, frozen = Path(root_s), Path(frozen_s)
    f = classifier(str(frozen))
    I = tuple(job['I']); m, n, index = job['m'], job['n'], job['index']
    counts = dict(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    unknown, raw_digest = [], hashlib.sha256()
    out = root / 'run' / 't54' / f'm{m:02}n{n:02}' / f'i{index:06}.jsonl.gz'
    if out.exists():
        raise RuntimeError(f'refusing to overwrite existing chunk {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_name(out.name + f'.{os.getpid()}.partial')
    with temp.open('wb') as rawfile:
        with gzip.GzipFile(fileobj=rawfile, mode='wb', compresslevel=1, mtime=0) as gz:
            for J in f.supports(n, CAP):
                counts['raw'] += 1
                if not f.canonical(I, J):
                    continue
                status, witness = f.classify(I, J)
                line = (json.dumps([I, J, status, witness], separators=(',', ':')) + '\n').encode()
                gz.write(line); raw_digest.update(line)
                counts['canonical'] += 1; counts[status] += 1
                if status == 'unknown':
                    unknown.append([I, J])
        rawfile.flush(); os.fsync(rawfile.fileno())
    os.replace(temp, out)
    event = {'event': 'chunk', 'total': 54, 'm': m, 'n': n, 'index': index, 'I': list(I),
             'counts': counts, 'unknown': unknown, 'raw_sha256': raw_digest.hexdigest(), 'gzip_sha256': sha(out)}
    receipt = {'schema': 1, 'job': {'total': 54, 'm': m, 'n': n, 'index': index, 'I': list(I)},
               'event_sha256': event_sha(event), 'gzip_sha256': event['gzip_sha256'],
               'input_manifest_sha256': manifest_sha, 'status': 'persisted_before_event'}
    atomic_json(root / 'receipts' / 't54' / f'm{m:02}n{n:02}' / f'i{index:06}.json', receipt)
    return event


def prepare(root: Path) -> None:
    if not Path('/Volumes/data').is_mount():
        raise RuntimeError('/Volumes/data is not mounted')
    if root.exists():
        raise RuntimeError(f'output already exists: {root}')
    for path, want in ((SRC / 'support54-plan.json', PLAN_SHA), (OLD / 'events.jsonl', EVENTS_SHA),
                       (OLD / 'received-subset-audit-gpu_partial39.json', AUDIT_SHA),
                       (SRC / 'frozen_support.py', FROZEN_SHA)):
        if sha(path) != want:
            raise RuntimeError(f'fixed input hash mismatch: {path}')
    plan = json.loads((SRC / 'support54-plan.json').read_text())
    if plan.get('format') != 'a30-total53-54-cap6-v1' or plan.get('job_total') != 185264:
        raise RuntimeError('unexpected full-plan contract')
    full = independent_full54(plan)
    received = load_received(full)
    missing = [job for key, job in full.items() if key not in received]
    if len(missing) != 90799:
        raise RuntimeError(f'expected 90799 missing total54 jobs, got {len(missing)}')
    root.mkdir(parents=True)
    (root / 'input').mkdir()
    for src, name in ((SRC / 'support54-plan.json', 'support54-plan.json'), (SRC / 'frozen_support.py', 'frozen_support.py'),
                      (OLD / 'events.jsonl', 'old-events.jsonl'), (OLD / 'received-subset-audit-gpu_partial39.json', 'received-audit.json'),
                      (AUDITOR, 'audited_chunk_verifier.py')):
        shutil.copy2(src, root / 'input' / name)
    manifest = {'schema': 1, 'scope': 'Only the missing total54 complement; not a finite exclusion conclusion.',
                'cap': CAP, 'total': 54, 'full_total54_jobs': len(full), 'historically_persisted_total54': len(received),
                'missing_total54_jobs': len(missing), 'source': {'full_plan_sha256': PLAN_SHA, 'old_events_sha256': EVENTS_SHA,
                'received_audit_sha256': AUDIT_SHA, 'frozen_support_sha256': FROZEN_SHA,
                'audited_verifier_sha256': sha(AUDITOR)}, 'missing_jobs': missing}
    atomic_json(root / 'input-manifest.json', manifest)
    append_event(root / 'run' / 'events.jsonl', {'event': 'local54_complement_started', 'input_manifest_sha256': sha(root / 'input-manifest.json'),
                 'expected_missing': len(missing), 'historically_persisted': len(received), 'cap': CAP})
    print(json.dumps({'prepared': str(root), 'missing': len(missing), 'manifest_sha256': sha(root/'input-manifest.json')}, sort_keys=True))


def run(root: Path, workers: int, seconds: float, limit: int | None) -> None:
    manifest_path = root / 'input-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest['missing_total54_jobs'] != 90799 or manifest['cap'] != CAP:
        raise RuntimeError('manifest contract failure')
    if sha(root / 'input' / 'frozen_support.py') != FROZEN_SHA:
        raise RuntimeError('copied classifier changed')
    events_path = root / 'run' / 'events.jsonl'
    existing = set()
    for line_no, line in enumerate(events_path.read_text().splitlines(), 1):
        ev = json.loads(line)
        if ev.get('event') == 'chunk':
            key = (ev['total'], ev['m'], ev['n'], ev['index'], tuple(ev['I']))
            if key in existing: raise RuntimeError(f'duplicate output event at {line_no}')
            existing.add(key)
    jobs = []
    for job in manifest['missing_jobs']:
        key = (job['total'], job['m'], job['n'], job['index'], tuple(job['I']))
        if key not in existing: jobs.append(job)
    if len(existing) + len(jobs) != 90799:
        raise RuntimeError('output event set not a subset of manifest')
    started = time.monotonic(); submitted = finished = 0; failures = []
    append_event(events_path, {'event': 'local54_session_started', 'workers': workers, 'seconds_limit': seconds, 'limit': limit,
                 'existing': len(existing), 'pending': len(jobs), 'input_manifest_sha256': sha(manifest_path)})
    with futures.ProcessPoolExecutor(max_workers=workers) as pool:
        active = {}; cursor = 0
        while cursor < len(jobs) or active:
            timeup = time.monotonic() - started >= seconds or (limit is not None and submitted >= limit)
            while not timeup and cursor < len(jobs) and len(active) < 2 * workers:
                fut = pool.submit(produce_one, str(root), str(root / 'input' / 'frozen_support.py'), sha(manifest_path), jobs[cursor])
                active[fut] = jobs[cursor]; cursor += 1; submitted += 1
                timeup = time.monotonic() - started >= seconds or (limit is not None and submitted >= limit)
            if not active:
                break
            ready, _ = futures.wait(active, timeout=1, return_when=futures.FIRST_COMPLETED)
            for fut in ready:
                job = active.pop(fut)
                try:
                    ev = fut.result()
                    append_event(events_path, ev); finished += 1
                except Exception as exc:
                    failures.append({'job': job, 'error': repr(exc)})
            if failures:
                break
    total_events = sum(1 for line in events_path.read_text().splitlines() if json.loads(line).get('event') == 'chunk')
    status = 'complete' if total_events == 90799 and not failures else 'partial' if not failures else 'failed'
    receipt = {'schema': 1, 'status': status, 'workers': workers, 'seconds_limit': seconds, 'limit': limit,
               'submitted': submitted, 'finished': finished, 'persisted_chunks_total': total_events,
               'remaining': 90799 - total_events, 'failures': failures, 'workers_joined': True,
               'input_manifest_sha256': sha(manifest_path), 'events_sha256': sha(events_path), 'elapsed_seconds': time.monotonic() - started}
    atomic_json(root / 'run' / 'last-session-receipt.json', receipt)
    append_event(events_path, {'event': 'local54_session_' + status, **receipt})
    print(json.dumps({'status': status, 'submitted': submitted, 'finished': finished, 'persisted_total': total_events,
                      'remaining': 90799 - total_events, 'failures': len(failures)}, sort_keys=True), flush=True)
    if failures:
        raise SystemExit(2)


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prepare'); p.add_argument('--root', type=Path, default=OUT_DEFAULT)
    r = sub.add_parser('run'); r.add_argument('--root', type=Path, default=OUT_DEFAULT); r.add_argument('--workers', type=int, default=4); r.add_argument('--seconds', type=float, default=14400); r.add_argument('--limit', type=int)
    args = ap.parse_args()
    if args.cmd == 'prepare': prepare(args.root)
    else:
        if args.workers != 4 or not (args.seconds > 0): raise SystemExit('local production requires exactly four workers and positive duration')
        run(args.root, args.workers, args.seconds, args.limit)

if __name__ == '__main__': main()
