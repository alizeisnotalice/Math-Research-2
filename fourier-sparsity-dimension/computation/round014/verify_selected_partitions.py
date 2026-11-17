#!/usr/bin/env python3
"""Independent exact certificate check for selected R14 DP partitions."""
import argparse, hashlib, json
from pathlib import Path

FULL=(1<<64)-1
def rank(vectors):
    piv={}
    for a in vectors:
        while a:
            j=a.bit_length()-1
            if j in piv:a^=piv[j]
            else:piv[j]=a;break
    return len(piv)
def words(basis):
    out=[0]
    for g in basis:out += [x^g for x in out]
    return out

p=argparse.ArgumentParser();p.add_argument('selected',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
raw=a.selected.read_bytes(); data=json.loads(raw); checked=[]
for record in data['boundary_records']:
    task=record['task']; basis=task['lift_basis']; d=task['lift_dimension']
    if rank(basis)!=d:raise AssertionError(('dependent',record['id']))
    code=words(basis)
    if len(code)!=(1<<d) or FULL not in code:raise AssertionError(('bad code',record['id']))
    union=0; terms=[]
    for text in record['partition']:
        cell=int(text,16)
        if not cell or cell not in code or union&cell:raise AssertionError(('bad partition cell',record['id'],text))
        union|=cell; sub=[g for g in code if g&~cell==0]; subrank=rank(sub)
        if len(sub)!=(1<<subrank):raise AssertionError(('subcode',record['id'],text))
        terms.append(subrank+1)
    if union!=FULL or sum(terms)!=record['kappa']:raise AssertionError(('certificate value',record['id']))
    if record['kappa']>record['upper']:raise AssertionError(('counterexample in selected',record['id']))
    checked.append({'id':record['id'],'D':d,'kappa':record['kappa'],'upper':record['upper'],'partition_terms':terms})
a.out.write_text(json.dumps({'status':'passed','input_sha256':hashlib.sha256(raw).hexdigest(),'records':checked},indent=2)+'\n')
