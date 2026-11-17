#!/usr/bin/env python3
"""All unital codes on four points; arbitrary algebra/core plus Q2 selector and diagonal addback."""
import argparse,hashlib,itertools,json,subprocess,time
from pathlib import Path
def span(b):
    s={0}
    for a in b:s|={v^a for v in list(s)}
    return s
def canonical(vs):
    piv={}
    for x in vs:
        for j in sorted(piv,reverse=True):
            if x>>j&1:x^=piv[j]
        if x:
            j=x.bit_length()-1
            for k in piv:
                if piv[k]>>j&1:piv[k]^=x
            piv[j]=x
    return tuple(piv[j] for j in sorted(piv,reverse=True))
def partitions(n):
    def rec(a):
        if len(a)==n:yield a;return
        for v in range(1+max(a,default=-1)+1):yield from rec(a+[v])
    yield from rec([])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();start=time.monotonic()
    parts=[tuple(sum(1<<j for j,v in enumerate(p) if v==c) for c in range(max(p)+1)) for p in partitions(4)]
    assert len(parts)==15
    codes={tuple(sorted(span([15]+[v for v in range(1,8) if s>>(v-1)&1]))) for s in range(128)}
    assert len(codes)==16
    results=[];cache={}
    for ct in sorted(codes):
        C=set(ct);Cb=canonical(C)
        legal=[cells for cells in parts if all(c in C for c in cells)]
        kapC=max(sum((sum(g&~cell==0 for g in C).bit_length()-1)+1 for cell in cells) for cells in legal)
        for cells in legal:
            Bb=canonical([g for g in C if all((g&e) in C for e in cells)])
            basis=[]
            for w in range(2):
                for g in Bb:basis.append(sum(((g>>(p&3))&1)<<p for p in range(64) if (p>>3)&1==w))
            for e in cells:basis.append(sum((((e>>(p&3))&1)*((p>>2)&1))<<p for p in range(64)))
            for g in Cb:basis.append(sum(((g>>(p&3))&1)<<p for p in range(64)))
            basis=canonical(basis)
            assert len(basis)<=12
            if basis not in cache:
                r=subprocess.run(['nice','-n','10',str(a.binary.resolve()),str(len(basis)),*map(str,basis)],capture_output=True,text=True,check=True,timeout=30)
                dp=json.loads(r.stdout);cw=span(basis);union=0;score=0
                for hexc in dp['partition']:
                    c=int(hexc,16);assert c in cw and not c&union;union|=c
                    sub=[g for g in cw if g&~c==0];rb=canonical(sub)
                    assert len(sub)==1<<len(rb);score+=len(rb)+1
                assert union==(1<<64)-1 and score==dp['kappa']
                cache[basis]=dp
            dp=cache[basis]
            results.append({'C_basis':Cb,'C_dimension':len(Cb),'base_capacity':kapC,'E_cells':cells,'B_basis':Bb,'lift_basis_hex':[hex(x) for x in basis],'lift_dimension':len(basis),'lift_capacity':dp['kappa'],'predicted_upper':2*kapC,'inequality_passed':dp['kappa']<=2*kapC,'optimum_partition':dp['partition']})
    out={'status':'passed' if all(r['inequality_passed'] for r in results) else 'counterexample','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'base_codes':len(codes),'algebra_pairs':len(results),'distinct_lifts':len(cache),'elapsed_seconds':time.monotonic()-start,'records':results,'scope':'All four-point function codes, not necessarily bent. Finite evidence for diagonal-addback capacity inequality only.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','base_codes','algebra_pairs','distinct_lifts','elapsed_seconds']}))
if __name__=='__main__':main()
