#!/usr/bin/env python3
"""D4-combine r030 joint tables and cross-check exact r018 distance moments."""
from __future__ import annotations
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
FILES={0:4,1:8,2:4,6:4,7:4,12:1}

def main(ref_path, input_dir, output):
    joint=defaultdict(int); classes=[]
    for v,mult in FILES.items():
        x=json.loads((input_dir/f'gap-n5-v{v}.json').read_text())
        assert x['n']==5 and x['u128_checked'] and len(x['holes'])==1
        h=x['holes'][0]; assert h['v']==v
        A=B=E=0; dd=defaultdict(int)
        for L,d,e2,c in h['joint_H_L_d_e2']:
            c=int(c); joint[L,d,e2]+=mult*c
            w=L*c; A+=w;B+=c;E+=e2*w;dd[d]+=w
        assert A==int(h['A2_v_from_fibres']) and B==int(h['B_v_from_D']) and E==int(h['weighted_e2_total'])
        assert [dd[i] for i in range(9)]==list(map(int,h['weighted_distance_counts']))
        classes.append({'v':v,'xy':h['xy'],'orbit_size':mult,'A2_v':str(A),'B_v':str(B),'weighted_e2_v':str(E)})
    A=sum(L*c for (L,d,e2),c in joint.items());B=sum(joint.values());E=sum(L*e2*c for (L,d,e2),c in joint.items())
    dh=defaultdict(int)
    for (L,d,e2),c in joint.items():dh[d]+=L*c
    ref=json.loads(ref_path.read_text())
    assert A==int(ref['A2']) and E==int(ref['squared_euclidean_total'])
    assert {str(d):dh[d] for d in dh}=={k:int(v) for k,v in ref['distance_counts'].items()}
    byL=[]
    for L in range(1,25):
        raw=sum(c for (ell,d,e),c in joint.items() if ell==L)
        if not raw:continue
        a=L*raw;e=sum(L*e2*c for (ell,d,e2),c in joint.items() if ell==L)
        byL.append({'L':L,'H_raw':str(raw),'A2_weighted':str(a),'weighted_e2':str(e),
                    'A2_share':str(Fraction(a,A)),'e2_share':str(Fraction(e,E)),
                    'mean_e2_over_n2_given_L':str(Fraction(e,25*a))})
    def block(name, pred):
        a=sum(L*c for (L,d,e),c in joint.items() if pred(L));e=sum(L*e2*c for (L,d,e2),c in joint.items() if pred(L))
        return {'name':name,'A2_weighted':str(a),'weighted_e2':str(e),'A2_share':str(Fraction(a,A)),
                'e2_share':str(Fraction(e,E)),'mean_e2_over_n2':str(Fraction(e,25*a))}
    sparse=[[L,d,e2,str(c)] for (L,d,e2),c in sorted(joint.items())]
    out={'status':'passed','n':5,'scope':'exact finite D4-combined stable-fibre geometry; no asymptotic inference',
         'A2_total':str(A),'B_total':str(B),'weighted_e2_total':str(E),'Q_e2_over_n2':str(Fraction(E,25*A)),
         'weighted_distance_counts':{str(d):str(dh[d]) for d in sorted(dh)},'by_L':byL,
         'blocks':[block('short_L_le_5',lambda L:L<=5),block('long_L_ge_6',lambda L:L>=6)],
         'joint_H_L_d_e2':sparse,'classes':classes,
         'crosscheck':'A2, distance histogram, and squared-Euclidean total match r018 n5-d4-root-moments.json'}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--input-dir',type=Path,default=HERE);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.reference,a.input_dir,a.output)
