"""Independent integer Rabin verification of the H9/G5 supplied factors."""
import json, sys, hashlib
from pathlib import Path
from math import lcm
from finite_fields import _mul_mod, _gcd_mod, _irreducible_mod, _prime, req

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def stretch(a,p):
    c=[0]*((len(a)-1)*p+1)
    for i,x in enumerate(a):c[i*p]=x
    return c
def check_factors(rows, target, ell):
    req(_prime(ell),'prime')
    out=[1];ds=[]
    for r in rows:
        a=r['ascending']
        req(all(type(x) is int for x in a),'integer encoding')
        req(type(r['degree']) is int and r['degree']==len(a)-1,'degree')
        req(type(r['multiplicity']) is int and r['multiplicity']==1,'multiplicity')
        a=[x%ell for x in a]
        _irreducible_mod(a,ell);out=_mul_mod(out,a,ell);ds.append(len(a)-1)
    req(out==[x%ell for x in target],'factor product')
    derivative=[i*out[i]%ell for i in range(1,len(out))]
    req(_gcd_mod(out,derivative,ell)==[1],'squarefree')
    return sorted(ds)
def verify(data):
    g=[1,0,0,-1,0,1];p5=[1]*5;p15=[1,-1,0,1,-1,1,0,-1,1];p25=[int(i%5==0) for i in range(21)]
    req(data['G5_ascending']==g,'G5 target')
    req(mul(g,p5)==[1,1,1,0,0,0,0,0,1,1],'H9 product')
    req(mul(p5,p15)==stretch(p5,3),'Phi5 x3')
    req(p25==stretch(p5,5),'Phi5 x5')
    types={}
    for key,ell,expected in [('irreducible5',2,[5]),('type_1_4',5,[1,4]),('transposition_power',3,[2,3])]:
        rec=data[key];req(rec['ell']==ell,'specified Galois prime')
        ds=check_factors(rec['factors'],g,ell);req(ds==expected,'Galois cycle type');types[str(ell)]=ds
    parts={'G_x3':stretch(g,3),'G_x5':stretch(g,5),'Phi5':p5,'Phi15':p15,'Phi25':p25}
    cycles=[]
    for p,ell,M in [(3,307,4),(5,1201,6)]:
        table=data['G5_substitution_factors'][str(ell)]
        req(set(table)==set(parts),'all parts')
        deg={k:check_factors(table[k],v,ell) for k,v in parts.items()}
        target=[1]
        for v in parts.values():target=mul(target,v)
        tm=[x%ell for x in target];req(_gcd_mod(tm,[i*tm[i]%ell for i in range(1,len(tm))],ell)==[1],'whole target squarefree')
        req(ell%p==1 and deg['G_x'+str(p)].count(p)==2,'two own fiber cycles')
        req(all(d%p!=0 for k,ds in deg.items() for d in ds if k!='G_x'+str(p) or d!=p),'no other p-divisible cycles')
        others=[d for k,ds in deg.items() for d in ds if k!='G_x'+str(p) or d!=p]
        req(lcm(*others)==M and M%p!=0,'full residual lcm')
        cycles.append({'p':p,'ell':ell,'M':M,'degrees':deg})
    return {'verified':True,'G5_cycle_types':types,'localized_double_cycles':cycles}
if __name__=='__main__':
    source=Path(sys.argv[1]);result=verify(json.loads(source.read_text()))
    result['input_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
    result['verifier_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result,indent=2))
