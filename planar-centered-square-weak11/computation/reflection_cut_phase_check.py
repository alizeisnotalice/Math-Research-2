#!/usr/bin/env python3
"""Exact cut-phase condition for frozen R35/R41 reflection states."""
import argparse
import hashlib
import json
import platform
from fractions import Fraction as F
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--r035', type=Path, required=True)
p.add_argument('--r041', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a=p.parse_args()
expected={
    'r035':'5021c54704af45b3ef8289c67a245926473434640372a934b2f15c4485f0c026',
    'r041':'ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3',
}
records=[]
for tag,path in [('r035',a.r035),('r041',a.r041)]:
    raw=path.read_bytes(); digest=hashlib.sha256(raw).hexdigest()
    assert digest==expected[tag]
    d=json.loads(raw); px=F(d['Px']); pts=[[F(v) for v in q] for q in d['pts']]
    assert F(d['lam'])==1 and all(0<=q[0]<px for q in pts)
    xmax=max(q[0] for q in pts)
    margin=xmax+4*px-30
    assert margin>0
    records.append({'template':tag,'input_sha256':digest,'atoms':len(pts),
        'Px':str(px),'xmax':str(xmax),'xmax_plus_4Px_minus_30':str(margin),
        'margin_decimal_display':float(margin),'strict_positive':True,
        'N5_left_capture_cut_lower':str(margin)})
assert records[0]['Px']==records[1]['Px']
assert records[0]['xmax']==records[1]['xmax']
out={'status':'PASS','python':platform.python_version(),
     'arithmetic':'fractions.Fraction; exact strict comparisons',
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'records':records,
     'consequences':[
         'C>xmax implies C+4Px>30.',
         'c_N-30=(N-5)Px+(C+4Px-30)>(N-5)Px.',
         'Every source seen in the left seam window at r<15 lies in the last five columns.',
         'For N>=5 the seam source window lies strictly right of x=0, so the left cut is invisible.'
     ]}
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(out,indent=2)+'\n')
print('PASS',[(z['template'],z['xmax_plus_4Px_minus_30']) for z in records])
