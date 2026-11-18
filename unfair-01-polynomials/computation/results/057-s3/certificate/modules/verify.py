#!/usr/bin/env python3
"""Bind or freshly regenerate the complete degree-eight S={3} certificate."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
SOURCES=ROOT/'modules'
ROOT_SHA='6da30bcf5f93f0fc121923d1979a37d0508327a5401d10c4e8d86ceb41fd546a'
GROUP_SHA='5cf9a01d716e19c9bf5dac8bb43e90ce49c4d238a09b434ce89d09d85528ca36'
TREE_SHA='1952d4f37bc538467cbab5a486b6ce7dac313e9e550fb5dc4d1177dd957e638a'
MOTHERS=869396
STATES={'empty_forward':869395,'empty_reverse':1}
REVERSE_ONLY=[0,0,7,1,45,1,0,0,0]
COUNTS={'groups':1915327,'mothers':1523605945,'positive':864141765,'cayley':196098368,
        'outer_plus_reject':192995818,'outer_minus_reject':862290,'inner_plus_reject':918205,
        'inner_minus_reject':452659,'pass':MOTHERS}


def require(condition,message):
    if not condition:raise ValueError(message)


def shafile(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def js(path):
    return json.loads(path.read_bytes())


def manifest_verify():
    path=ROOT/'manifest.json'
    manifest=js(path)
    require(manifest['format']==1,'package format')
    seen=set()
    for row in manifest['files']:
        rel=Path(row['rel'])
        require(not rel.is_absolute() and '..' not in rel.parts and rel!=Path(path.name) and rel not in seen,'unsafe manifest path')
        actual=ROOT/rel
        require(actual.is_file() and actual.stat().st_size==row['bytes'] and shafile(actual)==row['sha256'],f'bytes of {rel}')
        seen.add(rel)
    check=list(ROOT.rglob('*'))
    require(not any(p.is_symlink() for p in check),'cert symlink')
    actual={p.relative_to(ROOT) for p in check if p.is_file() and p!=path}
    require(seen==actual,'manifest inventory')
    return shafile(path)


def chunks(index):
    top=hashlib.sha256()
    total=0
    for item in index['parts']:
        with gzip.open(ROOT/'data'/item['file'],'rb') as f:
            payload=f.read()
        content=payload.splitlines(keepends=True)
        require(len(content)==item['rows'] and hashlib.sha256(payload).hexdigest()==item['raw_sha256'],'part binding')
        top.update(payload)
        total+=len(content)
        yield from content
    require(total==index['rows'] and top.hexdigest()==index['raw_sha256'],'raw ledger binding')


def archive_verify():
    manifest=manifest_verify()
    files=js(ROOT/'data/index.json')
    roots, trees=files['roots'],files['trees']
    require((roots['raw_sha256'],trees['raw_sha256'])==(ROOT_SHA,TREE_SHA),'expected frozen ledger hashes')
    seen=Counter()
    layers=0
    previous=None
    for line, record in zip(chunks(roots),chunks(trees),strict=True):
        row=json.loads(line)
        root,status,steps,death,live,digest=json.loads(record)
        require((previous is None or previous<row) and row==root and len(row)==9 and len(digest)==64,'row identity, order or transcript digest')
        previous=row
        require(status in STATES,'unknown tree in archive')
        if status=='empty_reverse':
            require(row==REVERSE_ONLY and (steps,death,live)==(122,26,237),'reverse-only regression')
        else:
            require(steps==death and 1<=death<=96 and live>=1,'forward depth')
        seen[status]+=1
        layers+=steps
    require(dict(seen)==STATES and layers==8869759,'tree cover and layer count')
    original=js(ROOT/'original/dual-tree-summary.json')
    require(original['input_sha256']==ROOT_SHA and original['ledger_sha256']==TREE_SHA and
            original['source_sha256']==shafile(SOURCES/'s3_trees.py') and original['states']==STATES and
            original['mothers']==MOTHERS and original['layers']==layers and original['unresolved']==[],
            'original primary receipt')
    third=js(ROOT/'original/third-tree-full.json')
    require(third['reference_sha256']==TREE_SHA and third['states']==STATES and third['layers']==layers and
            third['mothers']==MOTHERS and third['reverse_only']==REVERSE_ONLY and len(third['fallbacks'])==6,
            'original third receipt')
    return {'manifest_sha256':manifest,'mothers':MOTHERS,'states':dict(seen),'layers':layers,
            'root_sha256':ROOT_SHA,'tree_sha256':TREE_SHA}


def sorted_root_output(raw,out):
    rows=[tuple(json.loads(line)) for line in raw.read_bytes().splitlines()]
    require(len(rows)==len(set(rows))==MOTHERS,'new root count and uniqueness')
    with out.open('xb') as f:
        for row in sorted(rows):
            f.write(json.dumps(row,separators=(',',':')).encode()+b'\n')


def fresh(options,archive):
    if options.out.exists() or options.out.is_symlink():
        raise SystemExit('refuse existing --out; choose a new directory')
    out=options.out.resolve()
    require(ROOT!=out and ROOT not in out.parents,'fresh output must be outside the certificate')
    out.mkdir(parents=True,exist_ok=False)
    logs=out/'logs'
    logs.mkdir()
    scratch=out/'scratch'
    scratch.mkdir()
    environment={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','TMPDIR':str(scratch),'TEMP':str(scratch),'TMP':str(scratch)}
    py=sys.executable
    start=time.monotonic()
    now=lambda: datetime.now(timezone.utc).isoformat()
    began=now()
    commands=[]

    def run(label,command,stdout=None):
        command=list(map(str,command))
        path_out=stdout or logs/(label+'.stdout.txt')
        path_err=logs/(label+'.stderr.txt')
        started=time.monotonic()
        print(f'{now()} {label}',flush=True)
        with path_out.open('xb') as a,path_err.open('xb') as b:
            with subprocess.Popen(command,cwd=out,env=environment,stdout=a,stderr=b) as proc:
                pid=proc.pid
                code=proc.wait()
        item={'label':label,'command':command,'pid':pid,'exit_code':code,'seconds':time.monotonic()-started,
              'stdout_sha256':shafile(path_out),'stderr_sha256':shafile(path_err)}
        (logs/(label+'.json')).write_text(json.dumps(item,sort_keys=True)+'\n')
        commands.append(item)
        require(code==0,f'{label} failed with exit code {code}')
        return path_out,path_err

    gmp=options.gmp_prefix.resolve()
    cache=options.clang_module_cache.resolve()
    require(cache.is_dir(),'set --clang-module-cache to an existing shared cache directory')
    cc=options.cxx
    compiler,_=run('compiler-version',[cc,'--version'])
    common=[cc,'-std=c++20','-O3','-fno-fast-math','-ffp-contract=off','-Wno-deprecated-literal-operator',
            '-fmodules-cache-path='+str(cache)]
    lib=[f'-I{gmp/"include"}',f'-L{gmp/"lib"}','-lgmpxx','-lgmp']
    screen=out/'screen'
    boundary=out/'boundary'
    budget=out/'budget'
    run('compile-screen',[*common,SOURCES/'s3_tail_screen.cpp',*lib,'-o',screen])
    run('compile-fixture',[*common,'-Wno-return-type',SOURCES/'exact_boundary.cpp',*lib,'-o',boundary])
    run('compile-budget',[*common,SOURCES/'s3_budget_tree.cpp','-o',budget])
    fixture,_=run('exact-boundary',[boundary])
    require(js(fixture)==js(ROOT/'original/exact_boundary.json'),'fixture outcome')
    groups=out/'groups.txt'
    run('groups',[py,SOURCES/'s3_tail_groups.py','--out',groups])
    require(shafile(groups)==GROUP_SHA,'group stream')
    raw=out/'raw-roots.jsonl'
    _,counts_log=run('root-screen',[screen,groups],stdout=raw)
    screen_counts=js(counts_log)
    require(all(screen_counts.get(key)==value for key,value in COUNTS.items()),'root coverage and first-code counts')
    roots=out/'roots.jsonl'
    sorted_root_output(raw,roots)
    require(shafile(roots)==ROOT_SHA,'fresh complete root survivors')
    trees=out/'dual'
    run('dual-trees',[py,SOURCES/'s3_trees.py','--input',roots,'--out',trees,
                      '--workers',options.workers,'--depth',96,'--state-cap',10000])
    require(shafile(trees/'rows.jsonl')==TREE_SHA,'dual tree ledger')
    native=out/'native-budget.raw'
    run('third-native',[budget,roots],stdout=native)
    third=out/'third-comparison.json'
    run('third-exact-compare',[py,SOURCES/'s3_budget_crosscheck.py','--native',native,
                              '--reference',trees/'rows.jsonl','--out',third])
    require(js(third)==js(ROOT/'original/third-tree-full.json'),'third tree plus arbitrary precision fallback')
    require(manifest_verify()==archive['manifest_sha256'],'certificate changed during replay')
    report={'status':'passed_fresh_full','architecture':platform.machine(),'python':sys.version,
            'compiler':compiler.read_text().splitlines()[0],'started':began,'finished':now(),
            'seconds':time.monotonic()-start,'workers':options.workers,'cached_cores':0,
            'groups_sha256':GROUP_SHA,'tree_summary_sha256':shafile(trees/'summary.json'),
            'third_summary_sha256':shafile(third),'binaries_sha256':{p.name:shafile(p) for p in (screen,boundary,budget)},
            'archive':archive,'screen_counts':screen_counts,'commands':commands}
    (out/'fresh-report.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('status','seconds','archive')},sort_keys=True))


def main():
    if not __debug__:raise SystemExit('refuse optimized Python: stage modules require assertions')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,help='fresh full output; omit for archive verification only')
    parser.add_argument('--workers',type=int,default=8,choices=range(1,14))
    parser.add_argument('--cxx',default='clang++')
    parser.add_argument('--gmp-prefix',type=Path,default=Path('/opt/homebrew'))
    parser.add_argument('--clang-module-cache',type=Path,default=Path('/Volumes/data/CodexTemp/cache/shared-clang-modules'))
    options=parser.parse_args()
    if options.out is not None and (options.out.exists() or options.out.is_symlink()):
        raise SystemExit('refuse existing --out; choose a new directory')
    archive=archive_verify()
    if options.out is None:
        print(json.dumps({'status':'passed_archive',**archive},sort_keys=True))
    else:
        fresh(options,archive)


if __name__=='__main__':
    main()
