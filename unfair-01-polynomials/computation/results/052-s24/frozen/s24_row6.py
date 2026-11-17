"""Necessary sixth-row filter invariant under integer c^2/f translations."""
import json,importlib.util,hashlib,time
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('pair_screen',ROOT/'s24_low.py');p=importlib.util.module_from_spec(s);s.loader.exec_module(p);m=p.m
def main(input_path=None,out_path=None):
    begin=time.monotonic();input_path=Path(input_path) if input_path is not None else ROOT/'s24-low.json';x=json.loads(input_path.read_text());rows=[];counts=Counter()
    for z in x['rows']:
        D=z['D'];a=m.omega(z['a_m'],z['a_n'],D);d=m.omega(z['d_m'],z['d_n'],D);g2,q2,g4,q4=[z[k] for k in ('g2','q2','g4','q4')]
        A2=m.add(a,(g2,0));B2=m.add(m.bar(a),(q2,0));A4=m.add(m.add(d,m.scale(g2,a)),(g4,0));B4=m.add(m.add(m.bar(d),m.scale(q2,m.bar(a))),(q4,0))
        cross=m.add(m.mul(A2,B4,D),m.mul(A4,B2,D));assert m.sg(cross,D)>0 and m.sg(m.add((1,0),m.neg(cross)),D)>0
        U=m.add(m.scale(g4,a),m.scale(g2,d));V=m.add(m.scale(q4,m.bar(a)),m.scale(q2,m.bar(d)));base=m.add(m.add(U,V),cross);assert base[1]==0 and base[0].denominator==1
        choices=[]
        for g6 in m.shifts(U,D):
            q6=1-int(base[0])-g6
            if m.unit(m.add(V,(q6,0)),D):choices.append((g6,q6))
        assert len(choices)<=1;counts['input']+=1
        if not choices:counts['rejected6']+=1;continue
        counts['survived6']+=1;g6,q6=choices[0];rows.append({**z,'rho6':g6,'sigma6':q6})
    out={'status':'complete_necessary_row6_filter','scope':'rho6=G6+f-c^2, sigma6=Q6+f-c^2. No c/e/f full mother coverage.','counts':dict(counts),'rows':rows,'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pair_source_sha256':hashlib.sha256((ROOT/'s24_low.py').read_bytes()).hexdigest(),'seconds':time.monotonic()-begin};Path(out_path if out_path is not None else ROOT/'s24-row6.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'}))
    return out
if __name__=='__main__':main()
