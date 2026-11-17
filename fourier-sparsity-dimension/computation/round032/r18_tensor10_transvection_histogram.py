#!/usr/bin/env python3
"""Exact transvection histogram for the explicit R18 n=6 graph seed tensored with u·v."""
import hashlib,json,platform,sys,time
from pathlib import Path

N=1<<10
R18=Path(__file__).resolve().parents[1]/'round018/graph-cases.json'

def dot(a,b): return (a&b).bit_count()&1

def fwht(a):
    a=list(a); h=1
    while h<N:
        for i in range(0,N,2*h):
            for j in range(i,i+h):
                x,y=a[j],a[j+h]
                a[j],a[j+h]=x+y,x-y
        h*=2
    return a

def main():
    t0=time.monotonic()
    source=json.loads(R18.read_text())
    case=next(c for c in source['cases'] if c['r']==0 and c['n']==6)
    f0=int(case['f0_hex'],16)
    # Coordinates: x in F_2^6, u in positions 6,7, v in positions 8,9.
    F=[]
    for z in range(N):
        x=z&63; u=(z>>6)&3; v=(z>>8)&3
        F.append(((f0>>x)&1)^dot(u,v))
    assert set(abs(x) for x in fwht([1-2*b for b in F]))=={32}
    hist={}; walsh_hist={}; candidates=[]; tested=0
    for direction in range(1,N):
        derivative=[F[x]^F[x^direction] for x in range(N)]
        W=fwht([1-2*b for b in derivative])
        assert W[0]==0
        for linear_form in range(1,N):
            if dot(linear_form,direction): continue
            tested+=1
            weight=(N+W[linear_form])//4
            assert 4*weight==N+W[linear_form]
            hist[weight]=hist.get(weight,0)+1
            walsh_hist[W[linear_form]]=walsh_hist.get(W[linear_form],0)+1
            if weight==64:
                candidates.append({'v':direction,'lambda':linear_form,
                                   'derivative_walsh_at_lambda':W[linear_form]})
    assert tested==522753 and sum(hist.values())==tested
    descriptor={'schema':'a28-r032-r18-graph-tensor10-input-v1',
                'r18_input_path':str(R18),'r18_input_sha256':hashlib.sha256(R18.read_bytes()).hexdigest(),
                'r18_case':{'r':0,'n':6,'f0_hex':case['f0_hex']},
                'tensor':'F(x,u,v)=F0(x)+u dot v, with x bits 0..5, u bits 6..7, v bits 8..9'}
    inp=Path(__file__).with_name('r18-tensor10-transvection-input.json')
    inp.write_text(json.dumps(descriptor,indent=2)+'\n')
    result={'schema':'a28-r032-r18-graph-tensor10-transvection-histogram-v1',
            'scope':'one explicit R18 graph seed F0 tensored with the two-bit quadratic u·v; all nonidentity invertible transvections only; finite-family result, no projection scan',
            'arithmetic':'exact GF(2) truth tables and integer FWHT',
            'weight_formula':'|{x:lambda(x)=1,D_vF(x)=1}|=(1024+W_{D_vF}(lambda))/4, since every nonzero bent derivative is balanced',
            'input_path':str(inp),'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),
            'transvections_tested':tested,
            'weight_histogram':{str(k):hist[k] for k in sorted(hist)},
            'derivative_walsh_histogram_nonzero_lambda':{str(k):walsh_hist[k] for k in sorted(walsh_hist)},
            'minimum_positive_weight':min(k for k in hist if k>0),
            'minimum_derivative_walsh_nonzero_lambda':min(walsh_hist),
            'maximum_derivative_walsh_nonzero_lambda':max(walsh_hist),
            'weight64_candidate_count':len(candidates),'candidates':candidates,
            'elapsed_seconds':time.monotonic()-t0,'python':sys.version,'platform':platform.platform(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output=Path(__file__).with_name('r18-tensor10-transvection-histogram.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['transvections_tested','weight_histogram','minimum_positive_weight','weight64_candidate_count','elapsed_seconds']},indent=2))

if __name__=='__main__': main()
