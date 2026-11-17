#!/usr/bin/env python3
"""Exact r030 geometry-table versus independent r029 L-tail cross-check."""
from __future__ import annotations
import hashlib,json
from pathlib import Path

HERE=Path(__file__).resolve().parent
FILES={0:'tail-full-n5-v0.json',1:'tail-full-n5-v1.json',2:'tail-full-n5-v2.json',
       6:'tail-full-n5-v6.json',7:'tail-full-n5-v7.json',12:'tail-full-n5-center.json'}
MULT={0:4,1:8,2:4,6:4,7:4,12:1}
def sha(p:Path)->str:return hashlib.sha256(p.read_bytes()).hexdigest()
def main(geo_path,r29,output):
    geo=json.loads(geo_path.read_text())
    cells={(L,d,e):int(c) for L,d,e,c in geo['joint_H_L_d_e2']}
    raw=[sum(c for (L,d,e),c in cells.items() if L==t) for t in range(1,25)]
    details=[]; inputs={}
    for v,name in FILES.items():
        p=r29/name; inputs[str(v)]={'path':str(p),'sha256':sha(p)}
        r=json.loads(p.read_text()); assert len(r['holes'])==1
        x=r['holes'][0]
        gp=geo_path.parent/f'gap-n5-v{v}.json'; g=json.loads(gp.read_text())['holes'][0]
        assert x['v']==v==g['v']
        assert x['histogram_L']==g['histogram_L'],(v,'hist')
        assert x['A2_v_from_fibres']==g['A2_v_from_fibres'],(v,'A')
        assert x['B_v_from_D']==g['B_v_from_D'],(v,'B')
        details.append({'v':v,'xy':g['xy'],'orbit_size':MULT[v],'complete_histogram_match':True,
                        'A2_v_match':True,'B_v_match':True})
    summary=r29/'n5-tail-summary.json'; inputs['r29_summary']={'path':str(summary),'sha256':sha(summary)}
    rsum=json.loads(summary.read_text())
    assert raw==list(map(int,rsum['L_histogram_global'])),('global L histogram')
    assert sum(raw)==int(rsum['B_total'])==int(geo['B_total'])
    assert sum((t+1)*z for t,z in enumerate(raw))==int(rsum['A2_total'])==int(geo['A2_total'])
    out={'status':'passed','scope':'six D4 puncture representatives and D4-combined full L histogram',
         'r29_inputs':inputs,'per_hole':details,'global_histogram_match':True,
         'global_B_match':True,'global_A2_match':True}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--geometry',type=Path,required=True);p.add_argument('--r29',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.geometry,a.r29,a.output)
