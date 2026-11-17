"""Independent root allocation audit using exact disks and fixed-point intervals.

No generator imports. The Taylor tail uses a different, deliberately coarse
bound D^3 (|Re c|+|Im c|+1)^D r^2 (binary C, 0<r<=1).
Coefficient propagation rounds outwards to a decimal integer grid. Only the
sign of the normalized coefficient is needed: its monic numerator and its
nonzero normalizing constant must have opposite signs.
"""
import argparse, hashlib, json, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
S=10**100

def unpack(z):
    assert int(z['d'])>0
    return F(int(z['n']),int(z['d']))

def cmul(a,b):
    return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]

def evaluate(coeff,z):
    v=(F(0),F(0))
    for a in reversed(coeff):
        w=cmul(v,z);v=(w[0]+a,w[1])
    return v

def outward(lo,hi):
    assert lo<=hi
    return lo.numerator*S//lo.denominator,-((-hi.numerator*S)//hi.denominator)

def imul(x,y):
    v=[a*b for a in x for b in y]
    return min(v)//S,-((-max(v))//S)

def convolution(p,q,degree):
    out=[(0,0)]*min(len(p)+len(q)-1,degree+1)
    for i,x in enumerate(p):
        for j,y in enumerate(q[:len(out)-i]):
            lo,hi=imul(x,y);u,v=out[i+j];out[i+j]=(u+lo,v+hi)
    return out

def group_factor(g,r):
    a,b=g[0]
    if len(g)==1:
        assert b==0
        return [outward(-a-r,-a+r),(S,S)]
    assert len(g)==2 and g[1]==(a,-b) and b!=0
    mid=a*a+b*b;rad=2*(abs(a)+abs(b))*r+r*r
    return [outward(mid-rad,mid+rad),outward(-2*a-2*r,-2*a+2*r),(S,S)]

def main():
    if not __debug__:
        raise RuntimeError('Run without -O: assertions are certificate gates')
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,required=True);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--diagnostics',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    start=time.monotonic();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    cert=json.loads(a.certificate.read_text());inputs=json.loads(a.input.read_text());diag=json.loads(a.diagnostics.read_text())
    assert cert['input_sha256']==sha(a.input) and cert['diagnostics_sha256']==sha(a.diagnostics)
    ids=list(range(18,24));assert sorted(cert['tasks'])==ids
    assert sorted(c['task'] for c in cert['certificates'])==ids
    taskmap={t['id']:t for t in inputs['tasks']};expected={t:[] for t in ids}
    for x in diag:
        if x['task'] in ids and x['classification']=='nearunknown':expected[x['task']].append(x['mask'])
    assert sum(map(len,expected.values()))==235
    for rows in expected.values():assert len(rows)==len(set(rows))
    r=unpack(cert['radius']);assert 0<r<=1
    checked=[];summary=[]
    for tc in cert['certificates']:
        tid=tc['task'];task=taskmap[tid];c=task['C'];D=len(c)-1
        assert D==54 and tc['degree']==54 and c[0]==c[-1]==1 and all(x in (0,1) for x in c)
        groups=[[(unpack(z['re']),unpack(z['im'])) for z in g] for g in tc['centres_by_group']]
        assert len(groups)==tc['groups']==task['groups']==len(task['roots'])
        for g,seeds in zip(groups,task['roots']):
            assert len(g)==len(seeds)
            for z,seed in zip(g,seeds):
                assert max(abs(z[i]-F(str(seed[i]))) for i in (0,1))<F(1,10**8)
        flat=[z for g in groups for z in g];assert len(flat)==D
        separation=min(max(abs(x[k]-y[k]) for k in (0,1)) for i,x in enumerate(flat) for y in flat[:i])
        assert separation>2*r
        margins=[]
        for z in flat:
            v=evaluate(c,z);dv=evaluate([i*c[i] for i in range(1,D+1)],z)
            tail=D**3*(abs(z[0])+abs(z[1])+1)**D*r*r
            margin=max(map(abs,dv))*r-sum(map(abs,v))-tail
            assert margin>0;margins.append(margin)
        factors=[group_factor(g,r) for g in groups]
        rows=tc['mask_certificates'];assert [q['mask'] for q in rows]==expected[tid]
        for row in rows:
            mask=row['mask'];assert type(mask) is int and 0<=mask<(1<<(len(groups)-1))
            bits=[0]+[(mask>>(len(groups)-1-i))&1 for i in range(1,len(groups))]
            assert bits==row['bits_by_group'] and row['negative_factor'] in ('A','B')
            side=1 if row['negative_factor']=='A' else 0;index=row['negative_index'];assert type(index) is int and index>0
            degree=sum(len(g) for b,g in zip(bits,groups) if b==side);assert index<=degree
            p=[(S,S)]
            for b,factor in zip(bits,factors):
                if b==side:p=convolution(p,factor,index)
            const=p[0];coef=p[index]
            assert (const[0]>0 and coef[1]<0) or (const[1]<0 and coef[0]>0)
            checked.append({'task':tid,'mask':mask,'side':row['negative_factor'],'index':index,'raw_constant_scaled':list(map(str,const)),'raw_coefficient_scaled':list(map(str,coef))})
        summary.append({'task':tid,'roots':D,'masks':len(rows),'min_coarse_rouche_margin':str(min(margins)),'min_center_separation':str(separation)})
    assert len(checked)==235
    result={'status':'passed_independent_coarse_rouche_and_outward_integer_intervals','scope':'Only 235 specified exact root allocations; no exclusion of other masks or general conjecture.','grid_denominator':str(S),'masks':235,'tasks':summary,'sign_witnesses':checked,'elapsed_seconds':time.monotonic()-start,'python':sys.version,'bindings':{k:sha(p) for k,p in [('certificate',a.certificate),('input',a.input),('diagnostics',a.diagnostics),('verifier',Path(__file__))]}}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['status','masks','elapsed_seconds']}))

if __name__=='__main__':main()
