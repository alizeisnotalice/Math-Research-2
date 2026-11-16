"""Fixed H8(2,3) exact certificate checker, independent of search/cache code.

Only integer/Fraction arithmetic. Rechecks root isolation, source mapping,
all canonical layouts, designated coefficient enclosures and standard class.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse, hashlib, json, math
from root_ball_common import (complex_from_json, complex_from_strings, parse_fraction,
    polynomial_rouche_bounds, lift_rouche_bounds, disk_disjoint, abs_bounds)
from verify_modular_galois import verify as verify_galois

H = [1,1,1,0,0,0,0,1,1]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def integer(x):
    assert type(x) is int
    return x

def convolve(a,b):
    return [sum(a[i]*b[k-i] for i in range(len(a)) if 0<=k-i<len(b))
            for k in range(len(a)+len(b)-1)]
def stretch(a,n):
    z=[0]*(n*(len(a)-1)+1); z[::n]=a; return z

def unique_by(rows,key):
    assert isinstance(rows,list)
    out={}
    for row in rows:
        k=key(row); assert k not in out, ('duplicate', k); out[k]=row
    return out

def read_disk(entry):
    d=entry['disk']; c=complex_from_json(d['center']); r=parse_fraction(d['radius'])
    assert r>0
    return c,r

def separate(disks):
    for i,(c,r) in enumerate(disks):
        assert abs(c[1])>r
        for z,s in disks[:i]: assert disk_disjoint(c,r,z,s)

def coefficient_interval(factors,ids,index):
    # Recompute only the needed low coefficients. Higher powers cannot feed
    # a lower coefficient in a polynomial product.
    lo=[F(1)]; hi=[F(1)]
    for group in sorted(ids):
        fac=factors[group]; n=min(index+1,len(lo)+2); nl=[]; nh=[]
        for k in range(n):
            lower=F(0); upper=F(0)
            for j in range(max(0,k-len(lo)+1), min(2,k)+1):
                i=k-j; a,b=lo[i],hi[i]; c,d=fac[j]
                products=(a*c,a*d,b*c,b*d)
                lower+=min(products); upper+=max(products)
            nl.append(lower); nh.append(upper)
        lo,hi=nl,nh
    assert index<len(lo)
    return lo[index],hi[index]

def verify(input_path,balls_path,certificate_path,galois_path):
    source=json.loads(Path(input_path).read_text()); balls=json.loads(Path(balls_path).read_text())
    cert=json.loads(Path(certificate_path).read_text()); galois=json.loads(Path(galois_path).read_text())
    assert len(source['candidates'])==1
    raw=source['candidates'][0]
    assert raw['H_coefficients_ascending']==H and (raw['a'],raw['b'],raw['degree'])==(2,3,40)
    assert H!=H[::-1]
    C=convolve(stretch(H,2),stretch(H,3))
    assert set(C)<={0,1} and sum(C)==25 and raw['C_coefficients_ascending']==C
    gr=galois['certificates'] if isinstance(galois,dict) else galois
    assert len(gr)==1 and gr[0]['coefficients_ascending']==H
    verified_group=verify_galois(gr[0])
    assert balls['H']==H and (balls['a'],balls['b'])==(2,3)
    assert balls['input_sha256']==sha(input_path)
    assert cert['format']=='H8-a2b3-negative-Q-compact-v2'
    assert cert['input_sha256']==sha(input_path) and cert['balls_sha256']==sha(balls_path)
    assert cert['H']==H and (cert['a'],cert['b'],cert['codes_total'])==(2,3,1482)
    scale=int(balls['scale']); assert scale>0
    rawgroups=unique_by(raw['groups'],lambda g:integer(g['index']))
    rawblocks=unique_by(raw['blocks'],lambda g:integer(g['block']))
    assert set(rawgroups)==set(range(20)) and set(rawblocks)==set(range(4))
    mothers=unique_by(balls['mothers'],lambda r:integer(r['block']))
    assert set(mothers)==set(range(4))
    mother_disks={}; all_mother_disks=[]
    for block,m in mothers.items():
        assert m['id']==f'mother-{block}' and m['kind']=='mother'
        z,r=read_disk(m); assert z[1]>r
        assert z==complex_from_strings(rawblocks[block]['parent_root_positive_imag'])
        lhs,rhs=polynomial_rouche_bounds(H,z,r,scale)
        assert lhs<rhs
        assert parse_fraction(m['lhs'])==lhs and parse_fraction(m['rhs'])==rhs
        mother_disks[block]=(z,r)
        all_mother_disks.extend([(z,r),((z[0],-z[1]),r)])
    separate(all_mother_disks)
    lifts=unique_by(balls['lifts'],lambda r:integer(r['source_group_index']))
    assert set(lifts)==set(rawgroups)
    identities=set(); lift_disks=[]; factors={}; by_family={}
    for idx,g in rawgroups.items():
        block=integer(g['block']); assert block in mothers
        family=g['family']; assert family in ('a_lift','b_lift')
        q=2 if family=='a_lift' else 3; slot=integer(g['lift_index']); assert 0<=slot<q
        key=(block,family,slot); assert key not in identities; identities.add(key)
        item=lifts[idx]
        assert (item['block'],item['family'],item['lift_index'],item['exponent'])==(block,family,slot,q)
        assert item['kind']=='lift' and item['id']==f'{family}-{block}-{slot}' and item['parent']==f'mother-{block}'
        z,r=read_disk(item)
        assert z==complex_from_strings(g['lift_root'])
        ac,ar=mother_disks[block]
        assert ac==complex_from_strings(g['parent_root_positive_imag'])
        lhs,rhs=lift_rouche_bounds(ac,ar,z,r,scale,q)
        assert lhs<rhs and parse_fraction(item['lhs'])==lhs and parse_fraction(item['rhs'])==rhs
        lift_disks.extend([(z,r),((z[0],-z[1]),r)])
        low,high=abs_bounds(z,scale); low-=r; high+=r; assert low>0
        factors[idx]=[(low*low,high*high),(-2*(z[0]+r),-2*(z[0]-r)),(F(1),F(1))]
        by_family[key]=idx
    assert len(lift_disks)==40
    separate(lift_disks)
    slots=[]
    for block in range(4):
        aa=[by_family[(block,'a_lift',j)] for j in range(2)]
        bb=[by_family[(block,'b_lift',j)] for j in range(3)]
        r=rawblocks[block]
        assert r['a_lift_groups']==aa and r['b_lift_groups']==bb and r['all_group_indices']==aa+bb
        slots.append(aa+bb)
    # All choices are independently generated here; no options cache is read.
    all_layouts={}
    local_counts={}
    for E in (2,3,4,5,6):
        if E==6:
            choice_sets=[[(0,1)],[(0,1),(2,3,4)],[(0,1),(2,3,4)],[(0,1),(2,3,4)]]
        else:
            uv=[(u,v) for u in range(3) for v in range(4) if 3*u+2*v==E]
            assert len(uv)==1
            u,v=uv[0]
            choices=[tuple(a)+tuple(2+i for i in b) for a in combinations(range(2),u) for b in combinations(range(3),v)]
            choices.sort(key=lambda selected:sum(1<<j for j in selected))
            choice_sets=[choices]*4
        radices=[len(c) for c in choice_sets]; total=math.prod(radices); local_counts[E]=total
        for code in range(total):
            z=code; digits=[]; selected=[]
            for block,choices in enumerate(choice_sets):
                digit=z%len(choices); z//=len(choices); digits.append(digit)
                selected.extend(slots[block][j] for j in choices[digit])
            assert z==0 and len(set(selected))==len(selected)
            all_layouts[(E,code)]=(digits,set(selected))
    assert local_counts=={2:81,3:16,4:81,5:1296,6:8} and len(all_layouts)==1482
    rec=unique_by(cert['negative_certificates'],lambda r:(integer(r['E']),integer(r['code'])))
    unknown=unique_by(cert['unknown'],lambda r:(integer(r['E']),integer(r['code'])))
    assert set(unknown)=={(6,0)} and not set(rec)&set(unknown)
    assert len(rec)==1481 and set(rec)|set(unknown)==set(all_layouts)
    standard=unknown[(6,0)]
    assert standard['digits']==[0,0,0,0] and standard['terminal_standard'] is True
    chosen=all_layouts[(6,0)][1]
    assert chosen=={i for i,g in rawgroups.items() if g['family']=='a_lift'}
    assert set(rawgroups)-chosen=={i for i,g in rawgroups.items() if g['family']=='b_lift'}
    # Rouché and the separate degree counts identify these exact monic products
    # as H(x^2), H(x^3), whose integer product C was checked above.
    assert len(chosen)*2==16 and (20-len(chosen))*2==24
    bits=integer(cert['denominator_power']); assert 1<=bits<=1024
    denominator=1<<bits; margins=[]
    for key,r in sorted(rec.items()):
        digits,selected=all_layouts[key]
        assert r['digits']==digits and r['side'] in ('A','B')
        ids=selected if r['side']=='A' else set(rawgroups)-selected
        k=integer(r['coefficient_index']); assert 1<=k<2*len(ids)
        left=F(int(r['lo_n']),denominator); right=F(int(r['hi_n']),denominator)
        assert left<=right<0
        lower,upper=coefficient_interval(factors,ids,k)
        assert left<=lower<=upper<=right, ('bad interval', key)
        margins.append(-right)
    return {'status':'all_1482_layouts_verified_exactly','root_disks':40,'nonreal_mother_disks':8,
            'negative_layouts':1481,'standard_layouts':1,'constant_factors':'handled analytically as empty/full root sets',
            'layout_counts':local_counts,'minimum_certified_negative_margin':str(min(margins)),
            'arithmetic':'integers and Fraction only; no endpoint pinning; no options cache',
            'Galois':verified_group,'input_sha256':sha(input_path),'balls_sha256':sha(balls_path),
            'certificate_sha256':sha(certificate_path),'galois_sha256':sha(galois_path),
            'verifier_sha256':sha(__file__),
            'mathematical_dependency':'018 root-product balance and 025 H8 root-layout coverage proof; standard class from exact root-family degree counts'}

def main():
    ap=argparse.ArgumentParser()
    for name in ('input','balls','certificate','galois','record'):ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args();out=verify(a.input,a.balls,a.certificate,a.galois)
    a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
