#!/usr/bin/env python3
"""MP100 numerical blocks/cache for M(x^a)M(x^b), parameterized even-degree M."""
from __future__ import annotations
import argparse, hashlib, itertools, json
from math import comb, gcd
from pathlib import Path
import mpmath as mp
import numpy as np

mp.mp.dps = 100

def sha(p: Path) -> str: return hashlib.sha256(p.read_bytes()).hexdigest()
def fmt(x: mp.mpf) -> str: return mp.nstr(x, 95)
def parse(s: str) -> list[int]:
    h=[int(x.strip()) for x in s.split(',')]
    if len(h)<3 or h[0]!=1 or h[-1]!=1 or (len(h)-1)%2: raise SystemExit('mother must be monic, constant-one, even degree')
    return h
def conv(a: list[int], b: list[int]) -> list[int]:
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): z[i+j]+=x*y
    return z
def stretch(h: list[int], q: int) -> list[int]:
    z=[0]*((len(h)-1)*q+1); z[::q]=h; return z
def convmp(a: list[mp.mpf], b: list[mp.mpf]) -> list[mp.mpf]:
    z=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): z[i+j]+=x*y
    return z
def quad(z: mp.mpc) -> list[mp.mpf]: return [abs(z)**2, -2*mp.re(z), mp.mpf(1)]
def packed(p: list[mp.mpf]) -> dict:
    d=np.asarray([float(x) for x in p], dtype=np.float64); hi=d.astype(np.float32)
    return {'hi':hi.tolist(),'lo':(d-hi.astype(np.float64)).astype(np.float32).tolist(),'float64':d.tolist()}
def resid(h: list[int], z: mp.mpc) -> mp.mpf:
    return abs(sum(c*z**i for i,c in enumerate(h)))/sum(abs(c)*abs(z)**i for i,c in enumerate(h))
def local(a:int,b:int,u:int,v:int):
    ans=[]
    for ca in itertools.combinations(range(a),u):
        for cb in itertools.combinations(range(b),v):
            ans.append((sum(1<<x for x in tuple(ca)+tuple(a+j for j in cb)),tuple(ca),tuple(cb)))
    return sorted(ans)
def prod(ids: list[int], groups: dict[int,dict]) -> list[mp.mpf]:
    p=[mp.mpf(1)]
    for i in ids: p=convmp(p,groups[i]['coefficients_mp'])
    return p
def franel(n:int,q:int) -> int: return sum(comb(q,k)**n for k in range(q+1))
def parents(h:list[int], source:Path|None):
    n=(len(h)-1)//2
    if source:
        c=json.loads(source.read_text())['candidates'][0]
        if c.get('H_coefficients_ascending')!=h: raise SystemExit('parent-root source mother mismatch')
        r=sorted((g for g in c['groups'] if g.get('type')=='original_quadratic'),key=lambda g:int(g['block']))
        if len(r)!=n or [int(g['block']) for g in r]!=list(range(n)): raise SystemExit('bad parent-root source blocks')
        return [mp.mpc(g['parent_root_positive_imag'][0],g['parent_root_positive_imag'][1]) for g in r], 'reused_parent_root_input', sha(source)
    r=mp.polyroots(list(reversed(h)),maxsteps=1000,error=False)
    p=sorted((z for z in r if mp.im(z)>0),key=lambda z:(float(mp.re(z)),float(mp.im(z))))
    if len(p)!=n: raise SystemExit('positive-imaginary mother root count mismatch')
    return p,'mpmath_polyroots',None

def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument('--a',type=int,required=True); ap.add_argument('--b',type=int,required=True)
    ap.add_argument('--mother',required=True); ap.add_argument('--mother-name',required=True); ap.add_argument('--parent-roots-input',type=Path)
    ap.add_argument('--mother-certificate',type=Path); ap.add_argument('--out-input',type=Path,required=True); ap.add_argument('--out-options',type=Path,required=True)
    x=ap.parse_args(); a,b=x.a,x.b
    if not(1<a<b and gcd(a,b)==1): raise SystemExit('require 1<a<b and gcd(a,b)=1')
    if not Path('/Volumes/data').is_mount(): raise SystemExit('/Volumes/data is not mounted')
    h=parse(x.mother); n=(len(h)-1)//2; roots,root_source,root_source_sha=parents(h,x.parent_roots_input)
    if any(mp.im(z)<=0 or resid(h,z)>=mp.mpf('1e-75') for z in roots): raise SystemExit('invalid numerical parent root')
    groups=[]; blocks=[]
    for block,alpha in enumerate(roots):
        ai=[]; bi=[]
        for family,count,saved in [('a_lift',a,ai),('b_lift',b,bi)]:
            radius=abs(alpha)**(mp.mpf(1)/count); theta=mp.arg(alpha)
            for lift in range(count):
                z=radius*mp.e**(1j*(theta+2*mp.pi*lift)/count); index=len(groups)
                groups.append({'index':index,'block':block,'family':family,'lift_index':lift,'coefficients_mp':quad(z),'coefficients_ascending':[float(q) for q in quad(z)],'lift_root':[fmt(mp.re(z)),fmt(mp.im(z))],'parent_root_positive_imag':[fmt(mp.re(alpha)),fmt(mp.im(alpha))]}); saved.append(index)
        blocks.append({'block':block,'parent_root_positive_imag':[fmt(mp.re(alpha)),fmt(mp.im(alpha))],'a_lift_groups':ai,'b_lift_groups':bi,'all_group_indices':ai+bi})
    C=conv(stretch(h,a),stretch(h,b))
    if set(C)-{0,1}: raise SystemExit('M(x^a)M(x^b) not binary')
    terminal=a*b; balance={}
    for u in range(a+1):
        for v in range(b+1):
            E=b*u+a*v
            if 0<E<terminal:
                if E in balance: raise SystemExit(f'nonunique balance E={E}')
                balance[E]=(u,v)
    es=sorted(balance); byid={g['index']:g for g in groups}; cached=[]
    for block in blocks:
        options={}
        for E in es:
            u,v=balance[E]; rows=[]
            for mask,ca,cb in local(a,b,u,v):
                selected=[block['a_lift_groups'][j] for j in ca]+[block['b_lift_groups'][j] for j in cb]
                other=[j for j in block['all_group_indices'] if j not in set(selected)]
                rows.append({'mode':'mixed','u':u,'v':v,'selected_local_mask_decimal':str(mask),'selected_group_indices':selected,'selected_local_pairs':[{'family':'a_lift','lift_index':j,'group_index':block['a_lift_groups'][j]} for j in ca]+[{'family':'b_lift','lift_index':j,'group_index':block['b_lift_groups'][j]} for j in cb],'Ablock':packed(prod(selected,byid)),'Bblock':packed(prod(other,byid))})
            assert len(rows)==comb(a,u)*comb(b,v); options[str(E)]=rows
        ai,bi=block['a_lift_groups'],block['b_lift_groups']
        options[str(terminal)]=[{'mode':'all_a','u':a,'v':0,'selected_group_indices':ai,'selected_local_pairs':[{'family':'a_lift','lift_index':j,'group_index':ai[j]} for j in range(a)],'Ablock':packed(prod(ai,byid)),'Bblock':packed(prod(bi,byid))},{'mode':'all_b','u':0,'v':b,'selected_group_indices':bi,'selected_local_pairs':[{'family':'b_lift','lift_index':j,'group_index':bi[j]} for j in range(b)],'Ablock':packed(prod(bi,byid)),'Bblock':packed(prod(ai,byid))}]
        cached.append({'block':block['block'],'a_lift_groups':ai,'b_lift_groups':bi,'all_group_indices':block['all_group_indices'],'options':options})
    portable=[{k:v for k,v in g.items() if k!='coefficients_mp'} for g in groups]
    internal=sum((comb(a,balance[E][0])*comb(b,balance[E][1]))**n for E in es); fa,fb=franel(n,a),franel(n,b)
    assert internal==fa*fb//2-2; total=internal+2**(n-1); assert total==fa*fb//2+2**(n-1)-2
    raw={'name':f'coprime_{x.mother_name}_a{a}_b{b}','kind':'coprime_substitution','H_coefficients_ascending':h,'mother_name':x.mother_name,'a':a,'b':b,'degree':len(C)-1,'C_coefficients_ascending':C,'C_sha256':hashlib.sha256(json.dumps(C,separators=(',',':')).encode()).hexdigest(),'root_precision_decimal_digits':mp.mp.dps,'mother_root_max_relative_residual':fmt(max(resid(h,z) for z in roots)),'parent_roots_source':root_source,'parent_roots_source_sha256':root_source_sha,'mother_certificate':str(x.mother_certificate) if x.mother_certificate else None,'mother_certificate_sha256':sha(x.mother_certificate) if x.mother_certificate else None,'parent_block_count':n,'total_groups':len(portable),'groups':portable,'blocks':blocks,'balance_terminal_E':terminal,'balance_equation':'E=b*u+a*v; all parent blocks have identical E','E_options':[{'E':E,'u':balance[E][0],'v':balance[E][1],'per_block_options':comb(a,balance[E][0])*comb(b,balance[E][1])} for E in es]+[{'E':terminal,'special_modes':['all_a','all_b'],'per_block_options':2,'first_block_fixed_mode':'all_a','canonical_terminal_modes':2**(n-1)}],'balanced_enumeration':{'constant_and_full_excluded':True,'complement_canonical':True,'Franel_power':n,'Franel_a':fa,'Franel_b':fb,'internal_total':internal,'terminal_canonical_modes':2**(n-1),'formula':'F_n(a)F_n(b)/2 + 2^(n-1) - 2','nonconstant_complement_representatives':total},'status':'numerical_prepared_not_algebraically_certified'}
    x.out_input.parent.mkdir(parents=True,exist_ok=True); x.out_input.write_text(json.dumps({'generator':'parameterized MP100 coprime substitution','candidates':[raw]},indent=2)+'\n')
    cache={'format':'coprime-balanced-ds-block-options-v1','candidate':raw['name'],'kind':raw['kind'],'a':a,'b':b,'n_blocks':n,'balance_terminal_E':terminal,'mp_dps':mp.mp.dps,'block_count':n,'blocks':cached,'input_sha256':sha(x.out_input),'source_sha256':sha(Path(__file__)),'mixed_radix_order':'block order; first block radix is least significant; terminal fixes its all_a option','claim_boundary':'Numerical root blocks only; root balancing is a search restriction, not an exact sign certificate.'}
    x.out_options.write_text(json.dumps(cache,indent=2)+'\n'); print(json.dumps({'candidate':raw['name'],'a':a,'b':b,'n_blocks':n,'E_options':[r['E'] for r in raw['E_options']],'total':total,'input_sha256':sha(x.out_input),'options_sha256':sha(x.out_options)},indent=2))
if __name__=='__main__': main()
