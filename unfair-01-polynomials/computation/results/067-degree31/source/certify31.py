"""Generate or check degree-31 unit identities against original convolution equations."""
from argparse import ArgumentParser
from collections import Counter
import json
from pathlib import Path
import sys
import time

import sympy as sp
from sympy.polys.orderings import grevlex

sys.dont_write_bytecode=True
from census31 import DEGREE, OLD, TAIL_BITS, descriptor, sha
from hybrid_29_30 import kernel_basis, transform
from modular_core_29_30 import equations
from prefix_29_30 import model, slots
from trace_unit import unit_witness
from triangular_29_30 import forward_eliminate
from verify_degrees_17_18 import H
from verify_degrees_25_26 import MODS, check_modulus, pack, truncated_product, unpack

SOURCES=('certify31.py','census31.py','hybrid_29_30.py','modular_core_29_30.py',
         'prefix_29_30.py','trace_unit.py','triangular_29_30.py','verify_degree_19.py',
         'verify_degrees_17_18.py','verify_degrees_21_22.py','verify_degrees_23_24.py',
         'verify_degrees_25_26.py','verify_degrees_29_30.py')


def source_digests():
    here=Path(__file__).resolve().parent
    return {n:sha((here if n in ('certify31.py','census31.py') else OLD)/n) for n in SOURCES}


def original(early,mod):
    field,variables,E=equations(DEGREE,descriptor(early),mod)
    assert sum(x*y for x,y in zip(TAIL_BITS,mod[2],strict=True))%mod[0]==0
    return field,variables,E[:-1]+[E[-1]-variables[-1]]


def verify(R,E,h):
    assert len(E)==len(h)
    assert sum(x*y for x,y in zip(E,h,strict=True))==R.one


def unit(early,mod):
    R,vs,E=original(early,mod)
    pivots=forward_eliminate(R,vs,E,True)
    remaining=[f for f,h in pivots]
    if not remaining:
        return None,{'reason':'all_equations_zero'}
    for f,h in pivots:
        if f.is_ground:
            data=[u/f.const() for u in h]
            verify(R,E,data)
            return [pack(x,mod[0]) for x in data],{'method':'isolated_constant'}
    result=unit_witness(remaining)
    if result is not None:
        data=transform(result,[h for f,h in pivots],R)
        verify(R,E,data)
        return [pack(x,mod[0]) for x in data],{'method':'buchberger_unit'}
    used={i for i in range(len(vs)-1) if any(any(g[i] for g in f) for f in remaining)}
    ordered=tuple(vs[i].as_expr() for i in sorted(used))
    check=sp.groebner([f.as_expr() for f in remaining],*ordered,order='grevlex',modulus=mod[0])
    if len(check)!=1 or check[0]!=1:
        return None,{'reason':'proper_ideal_after_exact_substitutions','variables':list(map(str,ordered)),
                     'basis':[pack(R.from_expr(g),mod[0]) for g in check]}
    basis=kernel_basis(R,remaining,ordered,grevlex)
    for f,h in basis:
        if f.is_ground:
            data=transform([u/f.const() for u in h],[v for g,v in pivots],R)
            verify(R,E,data)
            return [pack(x,mod[0]) for x in data],{'method':'extended_groebner_unit'}
    raise ArithmeticError('unit ideal established, but explicit identity missing')


def exact_obstruction(early):
    shape=descriptor(early)
    t=sp.Symbol('t')
    Hpoly=sp.Poly.from_list(H[::-1],t,domain=sp.ZZ)
    free=sp.symbols(f'u0:{len(slots(DEGREE,shape))}')
    eq,A,B=model(t,DEGREE,shape,free)
    def coefficient_list(expression):
        if free:
            parts=sp.Poly(expression,*free).terms()
        else:
            parts=[((),expression)]
        out=[]
        for monomial,coefficient in parts:
            rem=sp.Poly(coefficient,t,domain=sp.ZZ).rem(Hpoly)
            if rem:
                out.append([monomial,[int(rem.nth(i)) for i in range(6)]])
        return out
    return {'variable_indices':slots(DEGREE,shape),
            'early_zero_equations':[[row,coefficient_list(g)] for row,g in eq],
            'truncated_product_rows':[[row,coefficient_list(truncated_product(A,B,row))] for row in range(24,31)]}


def make_record(early):
    failed=[]
    for mod in MODS:
        p=mod[0]
        h,meta=unit(early,mod)
        if h is not None:
            return {'shape':early,'status':'unit','prime':p,'multipliers':h,**meta,'earlier_moduli':failed}
        failed.append({'prime':p,**meta})
    return {'shape':early,'status':'unresolved','moduli':failed,'exact_system':exact_obstruction(early)}


def check_record(row):
    if row['status']=='unit':
        mod=next(m for m in MODS if m[0]==row['prime'])
        R,_,E=original(row['shape'],mod)
        verify(R,E,[unpack(x,R) for x in row['multipliers']])
        return
    assert row['status']=='unresolved'
    assert row['exact_system']==json.loads(json.dumps(exact_obstruction(row['shape'])))


def append(handle,value):
    handle.write(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n')


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    pa=ArgumentParser(description=__doc__)
    pa.add_argument('--census',type=Path,required=True)
    pa.add_argument('--start',type=int,default=0)
    pa.add_argument('--stop',type=int)
    pa.add_argument('--verify',type=Path)
    pa.add_argument('--out',type=Path,required=True)
    args=pa.parse_args()
    source=source_digests()
    census=json.loads(args.census.read_bytes())
    assert census['status']=='complete' and census['counts']['all']==27648
    assert all(source[n]==h for n,h in census['sources_sha256'].items())
    assert census['early_word']==[1]*12 and census['late_word']==[1]*7
    all_rows=census['ledger']['candidates']
    count=len(all_rows)
    lo=args.start
    hi=args.stop if args.stop is not None else count
    assert 0<=lo<=hi<=count
    for mod in MODS:
        check_modulus(*mod)
    census_hash=sha(args.census)
    started=time.monotonic()
    header={'status':'started','start':lo,'stop':hi,'census_sha256':census_hash,'sources_sha256':source,
            'python':sys.version,'sympy':sp.__version__,'mode':'verify' if args.verify else 'generate'}
    if args.verify:
        header['witness_sha256']=sha(args.verify)
    previous=args.verify.open() if args.verify else None
    accepted=Counter()
    if previous:
        prior=json.loads(next(previous))
        assert all(prior[k]==header[k] for k in ('start','stop','census_sha256','sources_sha256'))
        assert prior['mode']=='generate'
    with args.out.open('x',buffering=1) as output:
        append(output,header)
        for k in range(lo,hi):
            early=all_rows[k]['shape']
            if previous:
                row=json.loads(next(previous))
                assert row['index']==k and row['shape']==early
                check_record(row)
            else:
                row=make_record(early)
                row['index']=k
            append(output,row)
            accepted[row['status']]+=1
            if row['status']=='unit':
                accepted['mod'+str(row['prime'])]+=1
                accepted[row['method']]+=1
            else:
                print('UNRESOLVED',k,early,flush=True)
            if (k-lo+1)%100==0:
                print(k+1,'/',hi,dict(accepted),flush=True)
        receipt={'status':'complete','start':lo,'stop':hi,'counts':dict(accepted),'seconds':time.monotonic()-started}
        if previous:
            prior_end=json.loads(next(previous))
            assert prior_end['status']=='complete' and prior_end['counts']==receipt['counts'] and not previous.read(1)
            previous.close()
        append(output,receipt)
    print(receipt,flush=True)


if __name__=='__main__':
    main()
