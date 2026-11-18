"""Finite-field projections with affine pre-elimination and optional proof multipliers."""
from argparse import ArgumentParser
from collections import Counter
from hashlib import sha256
from pathlib import Path
from itertools import product
import json,sys,time

sys.dont_write_bytecode=True
import sympy as sp
from sympy.polys import distributedmodules as dm
from sympy.polys.orderings import lex
from probe_27_28 import slots,model,required_positive,proof_exists_positive,scout
from verify_degrees_25_26 import MODS,BITS,make_ring,truncated_product,pack,unpack,restrict,check_modulus


def candidate_bits(j,early):
    a,b=required_positive(early)
    return list(product(*((1,) if proof_exists_positive(a,b,k) else (0,1) for k in range(j-7,j))))


def equations(j,shape,mod):
    prime,t,weight=mod
    R,*vs=make_ring(prime,len(slots(j,shape)))
    z=vs[-1]
    E,A,B=model(R(t),j,shape,vs[:-1])
    rhs=sum(w*truncated_product(A,B,k) for k,w in enumerate(weight,j-7))
    return R,vs,[R(s) for _,s in E]+[z-rhs]


def substitution(R,f,i,L,x):
    """Return f(L) and an exact quotient of [f(x)-f(L)] by x-L."""
    total=R.zero
    quotient=R.zero
    powers=[R.one]
    for powers_xyz,c in f.items():
        e=powers_xyz[i]
        base=R.from_dict({powers_xyz[:i]+(0,)+powers_xyz[i+1:]:c})
        while len(powers)<=e:
            powers.append(powers[-1]*L)
        total+=base*powers[e]
        if e:
            quotient+=base*sum((x**(e-1-k))*powers[k] for k in range(e))
    assert f-total==quotient*(x-L)
    return total,quotient


def forward_eliminate(R,vars,E,track):
    pairs=[(g,[R.one if i==k else R.zero for k in range(len(E))] if track else None) for i,g in enumerate(E) if g]
    while True:
        unit=next((x for x in pairs if x[0].is_ground),None)
        if unit is not None:
            return [unit]
        best=None
        for k,(g,d) in enumerate(pairs):
            if any(sum(m)>1 for m in g):
                continue
            for m,v in g.items():
                if any(m[:-1]):
                    q=next(i for i,x in enumerate(m) if x)
                    rank=(len(g),k,q)
                    if best is None or rank < best[0]:
                        best=(rank,v)
        if best is None:
            return pairs
        (_,k,i),coef=best
        g,d=pairs.pop(k)
        u=R.domain.one/coef
        x=vars[i]
        L=-(g*u-x)
        factor=[f*u for f in d] if track else None
        revised=[]
        for h,provenance in pairs:
            if all(not e[i] for e in h):
                revised.append((h,provenance))
                continue
            y,q=substitution(R,h,i,L,x)
            if y:
                next_prov=[b-q*f for b,f in zip(provenance,factor,strict=True)] if track else None
                revised.append((y,next_prov))
        pairs=revised


def prove(j,shape,mod,track,witness=None):
    prime=mod[0]
    R,vars,E=equations(j,shape,mod)
    pair=[]
    if witness is not None:
        pair=[(unpack(rec['only_z'],R),[unpack(p,R) for p in rec['multipliers']]) for rec in witness]
    else:
        basic=forward_eliminate(R,vars,E,track)
        symbols=tuple(v.as_expr() for v in vars)
        if not basic:
            pass
        elif all(all(not any(m[:-1]) for m in f) for f,q in basic):
            pair=basic
        elif track:
            inputs=[dm.sdm_from_vector([g.as_expr()],lex,R.domain,gens=symbols) for g,q in basic]
            basis,trans=dm.sdm_groebner(inputs,dm.sdm_nf_mora,lex,R.domain,extended=True)
            for g,change in zip(basis,trans,strict=True):
                if any(any(m[1:-1]) for m,c in g):continue
                poly=R.from_expr(dm.sdm_to_vector(g,symbols,R.domain,n=1)[0])
                coef=[R.from_expr(p) for p in dm.sdm_to_vector(change,symbols,R.domain,n=len(basic))]
                converted=[sum(c*prov[k] for c,(_,prov) in zip(coef,basic,strict=True)) for k in range(len(E))]
                pair.append((poly,converted))
        else:
            G=sp.groebner([f.as_expr() for f,q in basic],*symbols,order='lex',modulus=prime)
            pair=[(R.from_expr(p.as_expr()),None) for p in G.polys if all(not p.degree(t) for t in symbols[:-1])]
    univ=[]
    saved=[]
    for f,m in pair:
        assert f and all(not any(e[:-1]) for e in f)
        if track:
            assert len(m)==len(E) and f==sum(x*y for x,y in zip(m,E,strict=True))
            saved.append({'only_z':pack(f,prime),'multipliers':[pack(x,prime) for x in m]})
        univ.append(tuple((e[-1],int(c)%prime) for e,c in f.items()))
    return saved,univ


def normalize(s):
    a,b,c=s
    return (a,b,tuple(c))


def run(j,raw,track,witness,progress=False):
    assumed={normalize(x['shape']):x for x in witness} if witness is not None else None
    record=Counter()
    results=[]
    for k,original in enumerate(raw['late'],1):
        shape=normalize(original)
        words=candidate_bits(j,shape[0])
        cnt=len(words)
        used={}
        record['shapes']+=1
        record['words']+=cnt
        for mod in MODS:
            if not words:break
            p=mod[0]
            replay=None if assumed is None else assumed[shape]['proofs'][str(p)]
            proofs,univ=prove(j,shape,mod,track,replay)
            before=len(words)
            words=restrict(words,mod,univ)
            record[f'mod{p}_shapes']+=1
            record[f'mod{p}_rejections']+=before-len(words)
            record[f'mod{p}_relations']+=len(univ)
            if track:used[str(p)]=proofs
        if words:
            record['live_shapes']+=1
            record['live_words']+=len(words)
            if progress:print(j,'live',shape,words,flush=True)
        if track or words:
            results.append({'shape':shape,'proofs':used,'survived':words,'word_count':cnt})
        if progress and k%100==0:
            print(j,k,'of',len(raw['late']),'live',record['live_shapes'],flush=True)
    if assumed is not None:assert len(assumed)==len(results)
    return {'degree':j,'early_counts':raw['counts'],'modular_counts':dict(record),'items':results}


def main():
    if not __debug__:raise SystemExit('run with assertions')
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('--out',type=Path,required=True)
    cli.add_argument('--witness',type=Path)
    cli.add_argument('--scout',type=Path,help='fast unsealed mode reads a previously generated early list')
    cli.add_argument('--only',type=int,choices=(27,28),help='test a single candidate degree')
    args=cli.parse_args()
    assert args.witness is None or args.scout is None
    if args.out.exists():raise SystemExit('refuse existing output')
    start=time.monotonic()
    elapsed_sections={}
    past=None if args.witness is None else {x['degree']:x['items'] for x in json.loads(args.witness.read_bytes())['degrees']}
    prior=None if args.scout is None else {x['j']:x for x in json.loads(args.scout.read_bytes())}
    for p,t,w in MODS:check_modulus(p,t,w)
    reports=[]
    for j in (27,28) if args.only is None else (args.only,):
        t0=time.monotonic()
        base=scout(j) if prior is None else prior[j]
        print(j,'early',base['counts'],flush=True)
        elapsed_sections[str(j)+'_early']=time.monotonic()-t0
        t0=time.monotonic()
        result=run(j,base,prior is None,None if past is None else past[j],progress=True)
        elapsed_sections[str(j)+'_modular']=time.monotonic()-t0
        reports.append(result)
        print(j,result['modular_counts'],flush=True)
    base=Path(__file__).parent
    names=('verify_degrees_17_18.py','verify_degree_19.py','verify_degrees_21_22.py','verify_degrees_23_24.py',
           'verify_degrees_25_26.py','probe_27_28.py',Path(__file__).name)
    data={'status':'complete','mode':'scout' if prior is not None else 'generate' if past is None else 'identities_only',
          'degrees':reports,'source_sha256':{n:sha256((base/n).read_bytes()).hexdigest() for n in names},
          'python':sys.version,'sympy':sp.__version__,'elapsed_seconds':time.monotonic()-start,'sections':elapsed_sections}
    if past is not None:data['witness_receipt_sha256']=sha256(args.witness.read_bytes()).hexdigest()
    with args.out.open('x') as f:f.write(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')


if __name__=='__main__':main()
