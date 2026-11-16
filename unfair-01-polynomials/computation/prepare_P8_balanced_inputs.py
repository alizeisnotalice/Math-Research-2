#!/usr/bin/env python3
"""High-precision mother-root lifting inputs for balanced H(x)H(x^L) search.

No high-degree root solve and no GPU.  Groups are numerical preparation only.
"""
from __future__ import annotations
import hashlib, json
from math import comb
from pathlib import Path
import mpmath as mp
import sympy as S

HERE=Path(__file__).resolve().parent
OUT=HERE/'stretch_balanced_candidates_L9_L11_20261115.json'
mp.mp.dps=100
LVALUES=(9,)
H=[1,0,1,1,1,0,1,1,1]

def conv(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def fmt(z): return mp.nstr(z,80)

def mother_roots():
    # Descending coefficients for x^6+x^5+x^2+x+1.
    rs=mp.polyroots(list(reversed(H)),maxsteps=1000,error=False)
    pos=sorted([z for z in rs if mp.im(z)>0],key=lambda z:(float(mp.re(z)),float(mp.im(z))))
    assert len(pos)*2==len(H)-1
    res=max(abs(sum(a*z**i for i,a in enumerate(H)))/sum(a*abs(z)**i for i,a in enumerate(H)) for z in rs)
    return pos,res

def group_from_z(z):
    return [float(abs(z)**2),float(-2*mp.re(z)),1.0]

def make(L,parents,res):
    h=H[:] # ascending; all-real-root exclusion certified separately
    hxL=[0]*((len(h)-1)*L+1)
    for i,a in enumerate(h): hxL[i*L]=a
    c=conv(h,hxL)
    assert set(c)<={0,1} and c[0]==c[-1]==1
    groups=[]; blocks=[]
    for bi,a in enumerate(parents):
        gidx=[]
        oi=len(groups); groups.append({'index':oi,'block':bi,'type':'original_quadratic','weight':L,
                                       'coefficients_ascending':group_from_z(a),
                                       'parent_root_positive_imag':[fmt(mp.re(a)),fmt(mp.im(a))]})
        gidx.append(oi)
        theta=mp.arg(a); radius=abs(a)**(mp.mpf(1)/L)
        lifts=[]
        for j in range(L):
            z=radius*mp.e**(1j*(theta+2*mp.pi*j)/L)
            gi=len(groups); groups.append({'index':gi,'block':bi,'type':'lift_quadratic','lift_index':j,'weight':1,
                                           'coefficients_ascending':group_from_z(z),
                                           'lift_root_from_positive_imag_parent':[fmt(mp.re(z)),fmt(mp.im(z))]})
            gidx.append(gi); lifts.append(gi)
        blocks.append({'block':bi,'parent_root_positive_imag':[fmt(mp.re(a)),fmt(mp.im(a))],
                       'original_group':oi,'lift_groups':lifts,'all_group_indices':gidx})
    rebuilt=[1.0]
    for g in groups:
        rebuilt=np_poly_mul(rebuilt,g['coefficients_ascending'])
    recon=max(abs(float(x)-y) for x,y in zip(rebuilt,c))
    # Complement canonical enumeration: E=1..L-1 lifts only; E=L first block original.
    per_e=[comb(L,e) for e in range(1,L)]
    nonconstant=sum(v**len(parents) for v in per_e)+2**(len(parents)-1)
    digest=hashlib.sha256(json.dumps(c,separators=(',',':')).encode()).hexdigest()
    return {'name':f'stretch_P8_L{L}','H_coefficients_ascending':h,'L':L,'degree':len(c)-1,
            'C_coefficients_ascending':c,'C_sha256':digest,'root_precision_decimal_digits':mp.mp.dps,
            'mother_root_max_relative_residual':fmt(res),'block_order':'increasing (real,imag) among positive-imag mother roots',
            'parent_block_count':len(parents),'total_groups':len(groups),'groups':groups,'blocks':blocks,
            'balanced_enumeration':{'E_range':[1,L],'for_E_1_to_Lminus1':'select E lift groups in every block; no original group',
              'for_E_L':'each block chooses original group or all L lift groups; fix block 0 to original',
              'complement_canonical':True,'constant_and_full_excluded':True,
              'nonconstant_complement_representatives':nonconstant},
            'numerical_group_product_inf_residual':recon,
            'status':'numerical_prepared_not_algebraically_certified'}

def np_poly_mul(a,b):
    c=[0.0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c

def main():
    global OUT
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True)
    OUT=parser.parse_args().out
    assert Path("/Volumes/data").is_mount()
    parents,res=mother_roots()
    records=[make(L,parents,res) for L in LVALUES]
    OUT.write_text(json.dumps({'generator':'mpmath mother roots then analytic Lth lifts','candidates':records},ensure_ascii=False,indent=2)+'\n')
    for r in records: print(r['name'],r['degree'],r['total_groups'],r['balanced_enumeration']['nonconstant_complement_representatives'],r['numerical_group_product_inf_residual'])

if __name__=='__main__': main()
