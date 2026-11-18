"""Run independently reconstructible prime channels for the n=12 height-five T."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import argparse
import json
import platform
import re
import subprocess

ROOT = Path(__file__).resolve().parent
CXX = ROOT/'rect_protected.cpp'
BIN = ROOT/'rect_protected'
CERT = ROOT/'pratt-61bit.json'


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def save(file, value):
    temp = file.with_suffix('.new')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(file)


def run(job, directory, primes):
    w,h,edges,index = job
    stem = f'{w}x{h}-{edges}-p{index}'
    stdout = directory/f'{stem}.json'
    stderr = directory/f'{stem}.resource'
    budget = 3000 if h == 5 else 32
    argv = ['/usr/bin/time','-l','./rect_protected',str(w),str(h),edges,
            str(index),'--max-table-mib',str(budget),'--progress']
    receipt = dict(job=job, command=argv, started_utc=now())
    with stdout.open('x') as raw, stderr.open('x') as resource:
        with subprocess.Popen(argv,cwd=ROOT,stdout=raw,stderr=resource) as process:
            receipt['time_pid'] = process.pid
            save(directory/f'{stem}.started.json', receipt)
            status = process.wait()
    receipt.update(waited_reaped=True, returncode=status, ended_utc=now())
    receipt.update(stdout_sha256=digest(stdout), resource_sha256=digest(stderr))
    save(directory/f'{stem}.finished.json',receipt)
    if status != 0:
        raise RuntimeError(f'{stem} exited {status}')
    counted = json.loads(stdout.read_text())
    assert (counted['width'], counted['height'], counted['unsafe_edges'], counted['prime_index']) == job
    assert counted['prime'] == primes[index]
    assert counted['processed_by_layer'] == counted['states_by_layer']
    rss = re.search(r'^\s*(\d+)\s+maximum resident set size$',stderr.read_text(),re.M)
    assert rss is not None
    return {**receipt, 'rss_bytes':int(rss[1]), 'planned_array_bytes':counted['max_table_bytes']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('name')
    parser.add_argument('index', type=int, nargs='+', choices=range(5))
    parser.add_argument('--bottom', action='store_true')
    args = parser.parse_args()
    assert len(args.index) == len(set(args.index))
    directory = ROOT/args.name
    directory.mkdir()
    bound = json.loads(CERT.read_text())['roots']
    state = dict(name=args.name, runner_pid=__import__('os').getpid(), host_arch=platform.machine(),
                 python=platform.python_version(), clang=subprocess.run(['clang++','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0],
                 source_sha256=digest(CXX),binary_sha256=digest(BIN),certificate_sha256=digest(CERT),started_utc=now())
    save(directory/'started.json',state)
    top = [(12,5,'B',index) for index in args.index]
    bottom = [(6,7,'TR',index) for index in range(4)] if args.bottom else []
    with ThreadPoolExecutor(max_workers=len(top)+len(bottom)) as batch:
        jobs = list(batch.map(lambda job: run(job,directory,bound),top+bottom))
    state.update(jobs=jobs, ended_utc=now(), waited_all=True)
    save(directory/'finished.json',state)


if __name__ == '__main__':
    main()
