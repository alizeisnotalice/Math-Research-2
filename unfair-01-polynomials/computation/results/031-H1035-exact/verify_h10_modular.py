"""Integer-only Rabin and cycle tests for the H10 two-fiber argument."""
import json,sys,hashlib
from pathlib import Path
from math import lcm
from finite_fields import _mul_mod,_gcd_mod,_irreducible_mod,_prime,req
H=[1,1,1,0,0,0,0,0,0,1,1]
def stretch(h,p):
    z=[0]*((len(h)-1)*p+1)
    for i,v in enumerate(h):z[i*p]=v
    return z
def factors(rows,target,p,key):
    req(_prime(p),'prime');product=[1];degrees=[]
    for row in rows:
        a=row[key];req(all(type(c)is int for c in a),'integer coefficients')
        req(type(row['multiplicity'])is int and row['multiplicity']==1,'multiplicity')
        if 'degree' in row:req(type(row['degree'])is int and row['degree']==len(a)-1,'degree')
        a=[c%p for c in a];_irreducible_mod(a,p);degrees.append(len(a)-1);product=_mul_mod(product,a,p)
    req(product==[c%p for c in target],'factor product')
    req(_gcd_mod(product,[i*product[i]%p for i in range(1,len(product))],p)==[1],'squarefree family')
    return sorted(degrees)
def verify(mother,cycles):
    req(len(mother['certificates'])==1,'one mother certificate');m=mother['certificates'][0]
    req(m['coefficients_ascending']==H and m['degree']==10 and cycles['H']==H,'fixed H10 target')
    req(m['sha256']==hashlib.sha256(json.dumps(H,separators=(',',':')).encode()).hexdigest(),'mother hash')
    types={r['prime']:factors(r['factors'],H,r['prime'],'coefficients_ascending') for r in m['reductions']}
    req(types=={7:[10],73:[1,9],19:[2,3,5]},'S10 types')
    req(len(cycles['certificates'])==2,'two fiber certificates');ans=[]
    for rec,(p,ell,M) in zip(cycles['certificates'],[(3,223,40),(5,3121,24)]):
        req(rec['target_p']==p and rec['ell']==ell and rec['target_family']==p,'cycle identity')
        req(set(rec['families'])=={'3','5'},'both families')
        ds={q:factors(rec['families'][str(q)],stretch(H,q),ell,'coefficients_ascending_mod_ell') for q in (3,5)}
        whole=_mul_mod(stretch(H,3),stretch(H,5),ell)
        req(_gcd_mod(whole,[i*whole[i]%ell for i in range(1,len(whole))],ell)==[1],'whole target squarefree')
        req(ell%p==1 and ds[p].count(p)==2,'two p-cycles')
        others=[d for q,v in ds.items() for d in v if q!=p or d!=p]
        req(all(d%p!=0 for d in others),'no other p-divisible cycle')
        req(lcm(*others)==M and M%p!=0 and rec['residual_degree_lcm_M']==M,'full residual lcm')
        ans.append({'p':p,'ell':ell,'M':M,'factor_degrees':ds})
    return {'status':'H10_S10_and_two_localized_fiber_rotations_exact_verified','mother_types':types,'cycles':ans}
if __name__=='__main__':
    m,c=map(Path,sys.argv[1:3]);out=verify(json.loads(m.read_text()),json.loads(c.read_text()))
    out.update(mother_input_sha256=hashlib.sha256(m.read_bytes()).hexdigest(),cycles_input_sha256=hashlib.sha256(c.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps(out,indent=2))
