#!/usr/bin/env python3
"""Regenerate large traces outside Git and invoke the independent verifier."""
import argparse, json, os, subprocess
from pathlib import Path
src=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--out',type=Path,default=Path('/Volumes/data/CodexTemp/T2/a022-r017-rank-reproduction'))
a=p.parse_args()
assert Path('/Volumes/data').is_mount(), 'data volume must be mounted'
out=a.out.resolve()
assert out.is_relative_to('/Volumes/data/CodexTemp'), 'use the data temporary volume'
out.mkdir(parents=True,exist_ok=True)
cache=Path('/Volumes/data/CodexTemp/cache/swift-modules');cache.mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env.update(CLANG_MODULE_CACHE_PATH=str(cache),TMPDIR=str(out),TMP=str(out),TEMP=str(out))
binary=out/'rank_audit'
subprocess.run(['swiftc','-O','-module-cache-path',str(cache),str(src/'rank_audit_final_source.swift'),'-o',str(binary)],env=env,check=True)
for r in json.loads((src/'trace-inputs.json').read_text()):
    subprocess.run([str(binary),str(r['n']),str(r['tau']),str(r['particles']),r['seed_hex'],r['method'],str(out/r['file'])],env=env,check=True)
subprocess.run(['python3',str(src/'verify_rank_traces.py'),str(out)],env=env,check=True)
