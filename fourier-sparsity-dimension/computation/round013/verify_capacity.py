#!/usr/bin/env python3
"""Run the exact C++ DP and independently verify optimum witnesses and bent cases."""
import argparse,hashlib,json,math,platform,subprocess,sys,time
from pathlib import Path
FULL=(1<<64)-1
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rank_basis(vs):
    piv={}
    for a in vs:
        while a:
            j=a.bit_length()-1
            if j in piv:a^=piv[j]
            else:piv[j]=a;break
    return list(piv.values())
def codewords(basis):
    c=[0]
    for g in basis:c += [x^g for x in c]
    return c
def mask(fn):return sum(int(fn(p))<<p for p in range(64))
def fwht(a):
    a=a[:];h=1
    while h<len(a):
        for j in range(0,len(a),2*h):
            for i in range(j,j+h):a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
        h*=2
    return a
def small_cases():
    out=[]
    for k in [1,2,3]:
        q=1<<k
        b=[mask(lambda p,j=j:((p>>k)&(q-1))==j) for j in range(q)]
        b += [mask(lambda p,j=j:(p>>j)&1) for j in range(k)]
        seed=[bin((p&(q-1))&((p>>k)&(q-1))).count('1')%2 for p in range(q*q)]
        out.append((f'MM_k{k}',b,seed,2*q))
    b=[mask(lambda p,j=j:((p>>2)&3)==j) for j in range(4)]
    b += [mask(lambda p:(p>>0)&1),mask(lambda p:(p>>1)&1),mask(lambda p:((p>>1)&1)*((p>>2)&1))]
    seed=[((p&1)*((p>>2)&1))^(((p>>1)&1)*((p>>3)&1)) for p in range(16)]
    out.append(('triangular_Q4',b,seed,8))
    base=[mask(lambda p,j=j:((p>>2)&3)==j) for j in range(4)]
    base += [mask(lambda p:p&1),mask(lambda p:(p>>1)&1)]
    b=[g&mask(lambda p,w=w:((p>>5)&1)==w) for w in [0,1] for g in base]
    b += [mask(lambda p:(p>>4)&1)]
    seed=[((p&1)*((p>>2)&1))^(((p>>1)&1)*((p>>3)&1))^(((p>>4)&1)*((p>>5)&1)) for p in range(64)]
    out.append(('MM_k2_Q2',b,seed,16))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--period-input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    started=time.monotonic();history=json.loads(a.period_input.read_text());cases=small_cases()
    for i,row in enumerate(history['rows']):
        b=[FULL]+[mask(lambda p,j=j:(p>>j)&1) for j in range(6)]+row['relative_generators']
        cases.append((f'period_{i}',b,None,None))
    records=[]
    for name,b,seed,expected in cases:
        d=len(b);assert rank_basis(b) and len(rank_basis(b))==d
        cw=codewords(b);assert len(set(cw))==1<<d and FULL in cw
        args=[str(d),*[str(v) for v in b]]
        begin=time.monotonic()
        proc=subprocess.run(['nice','-n','10',str(a.binary.resolve()),*args],capture_output=True,text=True,timeout=90,check=True)
        elapsed=time.monotonic()-begin;res=json.loads(proc.stdout);assert res['D']==d
        cells=[int(s,16) for s in res['partition']];union=0;detail=[]
        for cell in cells:
            assert cell and cell in cw and not(union&cell);union|=cell
            supported=[g for g in cw if g&~cell==0];basis=rank_basis(supported)
            assert len(supported)==1<<len(basis)
            detail.append({'mask':hex(cell),'size':cell.bit_count(),'subcode_dimension':len(basis),'subcode_basis':[hex(g) for g in basis]})
        assert union==FULL and sum(t['subcode_dimension']+1 for t in detail)==res['kappa']
        if expected is not None:assert res['kappa']==expected
        bent=None
        if seed is not None:
            n=len(seed);q=math.isqrt(n);assert q*q==n
            for g in cw:
                transform=fwht([1-2*(seed[i]^((g>>i)&1)) for i in range(n)])
                assert all(abs(x)==q for x in transform)
            bent={'native_truth_table_size':n,'all_phases_checked':len(cw),'Walsh_magnitude':q,'domain64_repetition':64//n}
        record={'case':name,'D':d,'kappa':res['kappa'],'basis_masks':[hex(g) for g in b],'input_args_sha256':hashlib.sha256((' '.join(args)+'\n').encode()).hexdigest(),'DP_exit_code':proc.returncode,'DP_elapsed_seconds':elapsed,'partition':detail,'bent_verification':bent}
        records.append(record)
        print(json.dumps({'case':name,'kappa':res['kappa'],'seconds':elapsed}),flush=True)
    assert len(records)==25 and all(r['kappa']==16 for r in records if r['case'].startswith('period_'))
    # Invalid masks: two repeated basis vectors; a valid code missing full support.
    negative=[]
    for args in [['2',str(FULL),str(FULL)],['1','1']]:
        p=subprocess.run([str(a.binary.resolve()),*args],capture_output=True,text=True,timeout=10)
        assert p.returncode!=0;negative.append({'args':args,'exit_code':p.returncode,'stderr':p.stderr.strip()})
    out={'status':'passed','source_sha256':sha(__file__),'cpp_sha256':sha(Path(__file__).with_name('capacity_dp.cpp')),'binary_sha256':sha(a.binary),'period_input_sha256':sha(a.period_input),'python':sys.version,'platform':platform.platform(),'compiler':subprocess.check_output(['clang++','--version'],text=True).splitlines()[0],'records':records,'negative_controls':negative,'elapsed_seconds':time.monotonic()-started,'scope':'Exact capacities for listed codes only. Witness checks prove lower bounds; exhaustive DP and its separately proved recurrence prove upper bounds. General selector equality is analytic.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
