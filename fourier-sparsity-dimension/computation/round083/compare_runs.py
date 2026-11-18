#!/usr/bin/env python3
"""Verify all GPU type histograms and hashes against the exact CPU census."""
import json
from pathlib import Path


def load(name):
    return json.loads(Path(name).read_text())


def match(gpu_names,cpu):
    minis,macbooks=(load(x) for x in gpu_names)
    checks=[]
    for A,B in zip(minis['cases'],macbooks['cases']):
        if A['name']!=B['name']:
            raise AssertionError('different fixture order')
        name=A['name']
        for al,bl in zip(A['levels'],B['levels']):
            depth=al['depth']
            if al['sha256']!=bl['sha256'] or al['n_kernels']!=bl['n_kernels'] or depth!=bl['depth']:
                raise AssertionError(('different GPU output',name,depth))
            c=next(l for l in cpu[name]['levels'] if l['depth']==depth)
            exact=[x for x in c['exact_classes'] if x['active_differences']>=depth]
            if exact!=al['exact_classes'] or exact!=bl['exact_classes']:
                raise AssertionError(('different CPU-GPU exact histogram',name,depth))
            checks.append(dict(name=name,depth=depth,kernels=al['n_kernels'],histogram_entries=len(exact),sha256=al['sha256']))
    return checks


def main():
    refs={z['name']:z for f in ('matching-mini.json','matching-macbook.json','guarded-odd-mini.json')
          for z in load(f)['cases']}
    groups=[('matching-metal-mini.json','matching-metal-macbook.json'),
            ('guarded-odd-metal-mini.json','guarded-odd-metal-macbook.json')]
    records=[r for group in groups for r in match(group,refs)]
    print(json.dumps(dict(status='pass',same_signed_kernel_digest_and_exact_cut_multiset=records,
                          kernels_per_device=sum(r['kernels'] for r in records)),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
