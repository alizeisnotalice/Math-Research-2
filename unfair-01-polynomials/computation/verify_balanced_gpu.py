#!/usr/bin/env python3
"""Independent verifier for round019 balanced GPU checkpoints.

Does not import the searcher.  It verifies prepared algebraic input exactly
where integer data allow, reconstructs floating group products, reimplements
the mixed-radix decoder, and NumPy-rechecks every retained CPU candidate.
Results remain numerical evidence only.
"""
from __future__ import annotations
import hashlib,json
from itertools import product
from math import comb
from pathlib import Path
import numpy as np



def iconv(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def opts(block,L,E):
    o=int(block['original_group']); lifts=[int(q) for q in block['lift_groups']]
    if E==L:return [1<<o,sum(1<<j for j in lifts)]
    out=[]
    for bits in range(1<<len(lifts)):
        if bits.bit_count()==E:
            out.append(sum(1<<lifts[q] for q in range(len(lifts)) if bits>>q&1))
    return out

def decode(code,options):
    mask=0
    for oo in options:
        mask|=oo[code%len(oo)];code//=len(oo)
    return mask,code

def npcheck(C,fs,mask):
    A=np.array([1.]);B=np.array([1.]); da=0
    for j,f in enumerate(fs):
        if mask>>j&1:A=np.convolve(A,f);da+=len(f)-1
        else:B=np.convolve(B,f)
    A=A/A[0];B=B/B[0]; conv=np.convolve(A,B)
    return {'degree_A':int(da),'degree_B':int(len(C)-1-da),'residual':float(abs(conv-C).max()),
            'minimum':float(min(A.min(),B.min())),'maximum':float(max(A.max(),B.max())),
            'leading_error':float(abs(A[-1]-1)+abs(B[-1]-1)),
            'fractionality':float(max(np.minimum(abs(A),abs(A-1)).max(),np.minimum(abs(B),abs(B-1)).max()))}

def verify_candidate(raw,ck,input_path,source_path):
    L=int(raw['L']);H=list(raw['H_coefficients_ascending']);C=list(raw['C_coefficients_ascending']);
    lifted=[0]*(6*L+1)
    for i,a in enumerate(H): lifted[i*L]=a
    exact=iconv(H,lifted)
    assert exact==C and set(C)<={0,1}
    assert hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest()==raw['C_sha256']
    fs=[np.asarray(g['coefficients_ascending'],dtype=float) for g in raw['groups']]
    product_poly=np.array([1.])
    for f in fs:product_poly=np.convolve(product_poly,f)
    group_err=float(abs(product_poly-np.asarray(C,float)).max()/max(1.,max(map(abs,C))))
    assert len(fs)==raw['total_groups']==3*(L+1)
    assert len(raw['blocks'])==3
    assert sorted(q for b in raw['blocks'] for q in [b['original_group']]+b['lift_groups'])==list(range(len(fs)))
    for bi,b in enumerate(raw['blocks']):
        oi=int(b['original_group']);li=[int(q) for q in b['lift_groups']]
        assert len(li)==L and raw['groups'][oi]['block']==bi and raw['groups'][oi]['weight']==L
        assert all(raw['groups'][q]['block']==bi and raw['groups'][q]['weight']==1 for q in li)
    expected=sum(comb(L,E)**3 for E in range(1,L))+4
    assert expected==raw['balanced_enumeration']['nonconstant_complement_representatives']
    assert raw['balanced_enumeration']['E_range']==[1,L]
    assert raw['balanced_enumeration']['complement_canonical'] is True
    assert raw['balanced_enumeration']['constant_and_full_excluded'] is True
    inp_sha=hashlib.sha256(input_path.read_bytes()).hexdigest();src_sha=hashlib.sha256(source_path.read_bytes()).hexdigest()
    coverage=[]; checked=0; max_field_delta=0.; max_recheck_res=0.; numerical_survivors=0
    for entry in ck['E_progress']:
        E=int(entry['E']); options=[opts(b,L,E) for b in raw['blocks']]
        if E==L:options[0]=[1<<int(raw['blocks'][0]['original_group'])]
        total=np.prod([len(q) for q in options],dtype=np.int64).item()
        cov=(entry['status']=='complete' and entry['total_codes']==total and entry['completed_code_half_open']==[0,total]
             and entry['counts']['codes_seen']==total)
        coverage.append({'E':E,'expected_codes':int(total),'covered':bool(cov),'closed_tolerance_fractional_masks':entry['counts']['closed_tolerance_fractional_masks']})
        for z in entry.get('top64',[]):
            mask,rest=decode(int(z['code']),options);assert rest==0 and mask==int(z['mask'])
            # Balance and canonical convention independently from code decoding.
            for bi,b in enumerate(raw['blocks']):
                oi=int(b['original_group']); li=[int(q) for q in b['lift_groups']]
                so=(mask>>oi)&1; sl=sum((mask>>q)&1 for q in li)
                if E<L: assert so==0 and sl==E
                else:
                    assert (so,sl) in ((1,0),(0,L))
                    if bi==0:assert (so,sl)==(1,0)
            assert mask not in (0,(1<<len(fs))-1)
            got=npcheck(np.asarray(C,float),fs,mask); checked+=1;max_recheck_res=max(max_recheck_res,got['residual'])
            numerical_survivors+=int(got['minimum']>=-1e-6 and got['maximum']<=1+1e-6 and got['leading_error']<=1e-6 and got['residual']<=1e-6 and got['fractionality']>1e-5)
            cpu=z.get('cpu_float64_recheck')
            if cpu:
                for a,b in [('degree_A','degree_A'),('degree_B','degree_B')]:assert got[a]==cpu[b]
                max_field_delta=max(max_field_delta,abs(got['residual']-cpu['convolution_max_abs_residual']),abs(got['minimum']-cpu['minimum_actual_coefficient']),abs(got['maximum']-cpu['maximum_actual_coefficient']),abs(got['leading_error']-cpu['normalized_leading_error']))
    assert [z['E'] for z in coverage]==list(range(1,L+1))
    assert sum(z['expected_codes'] for z in coverage)==expected
    assert ck['status']=='complete' and all(z['covered'] for z in coverage)
    assert ck['input_sha256']==inp_sha and ck['source_sha256']==src_sha
    assert max_field_delta<1e-12
    return {'saved_sample_numerical_survivors':numerical_survivors,'name':raw['name'],'L':L,'exact_H_times_HxL':True,'groups':len(fs),'group_product_relative_error_independent':group_err,
      'checkpoint_input_sha_matches':ck['input_sha256']==inp_sha,'checkpoint_source_sha_matches':ck['source_sha256']==src_sha,
      'expected_nonconstant_complement_representatives':expected,'coverage_sum':sum(z['expected_codes'] for z in coverage),
      'coverage':coverage,'all_E_complete':all(z['covered'] for z in coverage),'canonical_encoding_checked':'For E<L: no original and E lifts in each block; for E=L: original/all-lifts and block 0 original.',
      'retained_masks_rechecked':checked,'max_numpy_recheck_residual':max_recheck_res,'max_reported_cpu_field_delta':max_field_delta,
      'numerical_note':'No closed-tolerance fractional mask is a numerical search observation only, not a mathematical exclusion.'}

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('manifest',type=Path);args=ap.parse_args()
    base=args.manifest.parent;reports=[]
    for row in json.loads(args.manifest.read_text())['cases']:
        input_path=base/row['input'];source_path=Path(__file__).with_name('balanced_root_gpu.py')
        raw=json.loads(input_path.read_text())['candidates'][row['index']]
        ck=json.loads((base/row['checkpoint']).read_text())
        reports.append(verify_candidate(raw,ck,input_path,source_path))
    print(json.dumps({'scope':'completed L=9,10,11 balanced GPU checkpoints; L=12 not included','reports':reports},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
