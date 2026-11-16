"""Independent integer coverage/mask audit, including interrupted prefixes.

No searcher, cache, or numerical verifier is imported. Numeric saved-polynomial
reconstruction is a separate check; this script proves only ledger consistency.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(c, raw, source_hash):
    a,b=raw['a'],raw['b']; terminal=a*b
    assert math.gcd(a,b)==1 and 1<a<b
    assert c['source_sha256']==source_hash and c['candidate']==raw['name']
    H=raw['H_coefficients_ascending']; n=len(raw['blocks']); degree=len(H)-1
    assert all(type(x) is int and x in (0,1) for x in H) and H[0]==H[-1]==1
    assert n>=1 and degree==2*n
    S=[i for i,c in enumerate(H) if c]; target=[0]*(degree*(a+b)+1)
    for i in S:
        for j in S: target[a*i+b*j]+=1
    assert max(target)==1 and target==raw['C_coefficients_ascending']
    states={b*u+a*v:(u,v) for u in range(a+1) for v in range(b+1) if 0<b*u+a*v<terminal}
    keys=sorted(states)+[terminal]
    rows=c['E_progress']; observed=sorted(r['E'] for r in rows)
    assert len(set(observed))==len(observed) and observed==keys[:len(observed)]
    full=c['status']=='complete'
    assert not full or observed==keys
    assert len(raw['groups'])==n*(a+b)
    all_groups=[]
    for block in raw['blocks']:
        assert len(block['a_lift_groups'])==a and len(block['b_lift_groups'])==b
        all_groups.extend(block['a_lift_groups']+block['b_lift_groups'])
    assert sorted(all_groups)==list(range(n*(a+b)))
    expected_total=2**(n-1)-2+sum(math.comb(a,u)**n for u in range(a+1))*sum(math.comb(b,v)**n for v in range(b+1))//2
    count=samples=missing=0
    for entry in sorted(rows,key=lambda e:e['E']):
        E=entry['E']; options=[]
        if E<terminal:
            u,v=states[E]
            local=sorted([(tuple(aa)+tuple(a+j for j in bb)) for aa in itertools.combinations(range(a),u) for bb in itertools.combinations(range(b),v)],key=lambda x:sum(1<<j for j in x))
            radices=[len(local)]*n
            for block in raw['blocks']:
                groups=block['a_lift_groups']+block['b_lift_groups']
                options.append([[groups[j] for j in choice] for choice in local])
        else:
            radices=[1]+[2]*(n-1)
            options=[[block['a_lift_groups'],block['b_lift_groups']] for block in raw['blocks']]
            options[0]=options[0][:1]
        total=math.prod(radices); hi=entry['coverage_completed_half_open'][1]
        assert type(hi) is int and entry['coverage_completed_half_open']==[0,hi] and 0<=hi<=total
        assert entry['radices_block0_first']==radices and entry['total_codes']==total
        assert entry['counts']['codes_seen']==hi
        f=entry['counts']['fractional_masks']; inside=entry['counts']['closed_unit_interval_fractional_masks']
        assert type(f) is int and type(inside) is int and 0<=inside<=f<=hi
        if full or E!=observed[-1]: assert entry['status']=='complete' and hi==total
        if entry['status']=='complete': assert hi==total
        retained=entry['top64']; assert len(retained)==min(64,hi)
        seen=set()
        for row in retained:
            code=row['code']; assert type(code) is int and 0<=code<hi and code not in seen;seen.add(code)
            z=code;digits=[];selected=[]
            for radix,choices in zip(radices,options):
                digit=z%radix;z//=radix;digits.append(digit);selected.extend(choices[digit])
            assert z==0 and digits==row['mixed_radix_digits_block0_first']
            mask=sum(1<<i for i in selected)
            assert str(mask)==row['A_root_group_bitmask_decimal'] and hex(mask)==row['A_root_group_bitmask_hex']
            assert row['degree_A']==2*len(selected) and row['degree_B']==2*(len(raw['groups'])-len(selected))
            missing+=int('cpu_float64_independent_root_reconstruction' not in row)
        samples+=len(retained);count+=hi
    assert count<=expected_total and (not full or count==expected_total)
    return dict(a=a,b=b,blocks=n,full=full,codes=count,full_codes=expected_total,saved_rows=samples,missing_cpu_copies=missing)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive-dir',type=Path,required=True);ap.add_argument('--record',type=Path,required=True);ap.add_argument('--source-name',default='stream_coprime_ds_gpu_v2.py');args=ap.parse_args()
    root=args.archive_dir; mapping=json.loads((root/'archive-map.json').read_text())
    source_hash=sha(root/'scripts'/args.source_name); results=[]
    runs=mapping.get('runs')
    if runs is None:
        assert isinstance(mapping.get('files'),dict)
        runs=[]
        for rel in sorted(mapping['files']):
            if rel.startswith('checkpoints/'):
                tag=Path(rel).stem.removesuffix('-checkpoint')
                inp=f'inputs/{tag}-input.json'
                assert inp in mapping['files']
                runs.append({'archived_checkpoint':rel,'archived_input':inp})
        assert runs
    for run in runs:
        cp=root/run['archived_checkpoint']; ip=root/run['archived_input']
        c=json.loads(cp.read_text());raw=json.loads(ip.read_text())['candidates'][c['candidate_index']]
        assert sha(ip)==c['roots_input_sha256']
        result=audit(c,raw,source_hash); result.update(checkpoint_sha256=sha(cp),input_sha256=sha(ip));results.append(result)
    out=dict(status='integer_coverage_and_every_saved_mask_verified',runs=results,total_codes=sum(r['codes'] for r in results),total_saved_rows=sum(r['saved_rows'] for r in results),source_sha256=source_hash,auditor_sha256=sha(Path(__file__)),claim='Ledger and integer masks only; numerical arithmetic and unsaved coefficients are not certified.')
    args.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='runs'}))


if __name__=='__main__': main()
