"""Full O_K S24 necessary pair screen, using independent Fraction surds.

The root/Newton bounds give |a|,|abar|<55/8, |d|,|dbar|<419/8,
|2a^2-4d| in both embeddings <787/16. Only the fourth low row is
checked here: no claim of cofactor completion or mother-tail coverage.
"""
import importlib.util,json,hashlib,math,sys,time
from pathlib import Path
from collections import Counter
from fractions import Fraction as F
T=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials')
BASE=T/'sync/unfair-01-polynomials/computation/results/051-s12/frozen/independent_low.py'
s=importlib.util.spec_from_file_location('exact',BASE);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def main():
    assert __debug__;begin=time.monotonic();counts=Counter();rows=[];fields=[]
    for D in range(2,190):
        if any(D%(p*p)==0 for p in range(2,math.isqrt(D)+1)):continue
        aa=m.lattice(D,F(55,8),True)
        if not aa:continue
        dd=m.lattice(D,F(419,8),False);counts['fields']+=1;counts['a']+=len(aa);counts['d_per_field']+=len(dd);fields.append([D,len(aa),len(dd)])
        for am,an,a in aa:
            g2=-m.fl(a,D);A2=m.add(a,(g2,0));q2=int(1-2*a[0]-g2);B2=m.add(m.bar(a),(q2,0));assert m.unit(A2,D) and m.unit(B2,D) and m.add(A2,B2)==(1,0)
            for dm,dn,d in dd:
                counts['raw_pairs']+=1
                s4=m.add(m.scale(2,m.mul(a,a,D)),m.scale(-4,d))
                if not m.bound(s4,D,F(787,16)):counts['s4_rejected']+=1;continue
                counts['s4_pass']+=1
                A4base=m.add(d,m.scale(g2,a));B4base=m.add(m.bar(d),m.scale(q2,m.bar(a)))
                C4base=m.add(m.add(A4base,B4base),m.mul(A2,B2,D));assert C4base[1]==0 and C4base[0].denominator==1
                choices=[]
                for g4 in m.shifts(A4base,D):
                    q4=1-int(C4base[0])-g4;B4=m.add(B4base,(q4,0))
                    if m.unit(B4,D):choices.append((g4,q4))
                assert len(choices)<=1
                if not choices:counts['row4_rejected']+=1;continue
                counts['survived']+=1;g4,q4=choices[0];rows.append({'D':D,'a_m':am,'a_n':an,'d_m':dm,'d_n':dn,'g2':g2,'q2':q2,'g4':g4,'q4':q4})
    out={'status':'complete_necessary_S24_pair_and_row4_screen','scope':'All O_K pairs under s2/coefficient/s4 bounds and C2/C4 only; c,e,f and further rows not yet covered.','counts':dict(counts),'fields':fields,'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dependency_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest(),'python':sys.version,'seconds':time.monotonic()-begin}
    Path(__file__).with_name('s24-low.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('rows','fields')}))
if __name__=='__main__':main()
