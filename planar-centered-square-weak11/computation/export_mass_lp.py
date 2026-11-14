#!/usr/bin/env python3
"""Export exact conservative mass columns and two central-difference estimates."""
import argparse,json,os
from fractions import Fraction as F
from pathlib import Path
from periodic_halo import sha

def main():
    p=argparse.ArgumentParser();p.add_argument('--base',required=True,type=Path);p.add_argument('--target',required=True,type=Path);p.add_argument('--prior',required=True,type=Path);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    if not os.path.ismount('/Volumes/data'):raise RuntimeError('data not mounted')
    base=json.loads(args.base.read_text());target=json.loads(args.target.read_text());m=json.loads((args.prior/'manifest.json').read_text());rs={r['id']:r for r in json.loads((args.prior/'results.json').read_text())};n=len(base['pts']);columns={}
    if m['base_sha256']!=sha(args.base) or m['target_sha256']!=sha(args.target):raise ValueError('prior input mismatch')
    for task in m['tasks']:
        r=rs[task['id']]
        if 'error' in r:raise ValueError('prior task failure')
        for alias in task['aliases']:
            if alias.get('kind')=='baseline':continue
            key=(alias['quadrant'],tuple(alias['members']),alias['anchor']);row=columns.setdefault(key,{'aliases':[],'samples':{}})
            row['aliases'].append(alias);sk=(alias['eps'],alias['sign'])
            if sk in row['samples'] and row['samples'][sk]['id']!=r['id']:raise ValueError('duplicate direction conflict')
            row['samples'][sk]=r
    out=[]
    for (quad,members,anchor),row in sorted(columns.items()):
        vals=[F(p[2]) for p in target['pts']];ids=[quad*n+i for i in members];a=quad*n+anchor;G=sum(vals[i] for i in ids);q=min(G,vals[a]);vector=[F(0) for _ in vals]
        for i in ids:vector[i]=q*vals[i]/G
        vector[a]-=q
        if sum(vector)!=0:raise ValueError('column not conserved')
        estimates={}
        for eps in ('1/1000','1/10000'):
            plus=row['samples'][(eps,1)];minus=row['samples'][(eps,-1)];mp=(F(plus['inner_ratio'])+F(plus['outer_ratio']))/2;mm=(F(minus['inner_ratio'])+F(minus['outer_ratio']))/2;g=(mp-mm)/(2*F(eps));estimates[eps]={'gradient_float':float(g),'gradient_fraction_of_interval_midpoints':str(g),'plus_id':plus['id'],'minus_id':minus['id'],'plus_interval':[plus['inner_ratio'],plus['outer_ratio']],'minus_interval':[minus['inner_ratio'],minus['outer_ratio']]}
        out.append({'quadrant':quad,'members':list(members),'anchor':anchor,'unit_mass_increment':list(map(str,vector)),'estimates':estimates,'aliases':row['aliases']})
    lo=['0']*len(target['pts']);hi=[str(F(1001,1000)*F(base['pts'][i%n][2])) for i in range(len(target['pts']))]
    payload={'schema':'mass-finite-differences-lp-v1','scope':'central finite differences and numerical LP are candidate-generation heuristics, not derivatives or optimality certificates','base_input_path':str(args.base),'base_input_sha256':sha(args.base),'target_input_path':str(args.target),'target_input_sha256':sha(args.target),'source_manifest_sha256':sha(args.prior/'manifest.json'),'source_results_sha256':sha(args.prior/'results.json'),'export_script_sha256':sha(__file__),'base_label_masses':[str(F(base['pts'][i%n][2])) for i in range(len(target['pts']))],'target_masses':[str(F(p[2])) for p in target['pts']],'lower_mass_bound':lo,'upper_mass_bound':hi,'kappa_max':'1001/1000','budgets':[1,2,4,8,16,32,64,132],'trust_scales':['1/10000','3/10000','1/1000'],'columns':out};args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(payload,indent=2)+'\n');print('columns',len(out),'labels',len(lo),'output',args.output)
if __name__=='__main__':main()
