#!/usr/bin/env python3
"""Exact accounting and mixed-difference enclosures for the four R48 batches."""
import json,hashlib,sys
from fractions import Fraction as F
from pathlib import Path
from scan_supercell_2x2 import digest
from periodic_halo import sha
P=Path('/Volumes/data/CodexTemp/T2/planar-centered-square-weak11/round048/mass-geometry'); blobs={};allhash=set();records=0;new=0;reuse=0;summary=[]
for name in ['batch0','batch1','batch2','cross']:
 p=P/name;m=json.loads((p/'manifest.json').read_text());rs=json.loads((p/'results.json').read_text());s=json.loads((p/'summary.json').read_text());rr={r['id']:r for r in rs};assert len(rr)==len(rs)==len(m['tasks'])==s['completed'] and not s['errors'];assert s['results_sha256']==sha(p/'results.json');
 for t in m['tasks']:
  r=rr[t['id']];assert digest(t['input'])==t['input_sha256']==r['input_sha256']; assert F(r['tail_bound'])<1 and F(r['kappa'])>0 and F(r['delta'])>=0 and F(r['inner_ratio'])<=F(r['outer_ratio']);assert all(F(z[2])>0 for z in t['input']['pts']);allhash.add(t['input_sha256'])
 records+=s['points'];new+=s['new_calls'];reuse+=s['reused'];summary.append({'name':name,'points':s['points'],'new_calls':s['new_calls'],'reused':s['reused'],'errors':s['errors'],'manifest_sha256':sha(p/'manifest.json'),'results_sha256':sha(p/'results.json'),'summary_sha256':sha(p/'summary.json')});blobs[name]=(m,rr)
assert records==258 and new==204 and reuse==54 and len(allhash)==205
m,rr=blobs['cross'];lookup={(a['mass'],a['geometry']):rr[t['id']] for t in m['tasks'] for a in t['aliases']};mixed=[]
for g in ['q2_28_minus','q2_28_plus','all_26_minus','all_26_plus']:
 a,b,c,d=[lookup[x] for x in [('final_grid_winner',g),('final_grid_winner','baseline'),('baseline',g),('baseline','baseline')]]
 lo=F(a['inner_ratio'])-F(b['outer_ratio'])-F(c['outer_ratio'])+F(d['inner_ratio']);hi=F(a['outer_ratio'])-F(b['inner_ratio'])-F(c['inner_ratio'])+F(d['outer_ratio']);assert lo<=hi and (lo>0 or hi<0);mixed.append({'geometry':g,'lower':str(lo),'upper':str(hi),'sign':'positive' if lo>0 else 'negative'})
for mass in ['baseline','first_grid_winner','final_grid_winner']:
 b=lookup[(mass,'baseline')]
 for g in ['q2_28_minus','q2_28_plus','all_26_minus','all_26_plus']:assert F(lookup[(mass,g)]['outer_ratio'])<F(b['inner_ratio'])
out={'scope':'Exact bookkeeping and interval arithmetic; verifies stored certificates structurally, does not rerun halo kernels or establish a continuous family bound.','records':records,'unique_inputs':len(allhash),'new_calls':new,'reused_records':reuse,'batches':summary,'mixed_difference_definition':'(J(final mass,g)-J(final mass,baseline))-(J(original mass,g)-J(original mass,baseline))','mixed_differences':mixed,'all_12_geometry_changes_strictly_decrease_at_three_tested_masses':True,'script_sha256':sha(__file__),'python':sys.version};path=P/'exact_audit.json';path.write_text(json.dumps(out,indent=2)+'\n');print('audit',sha(path),'records',records,'unique',len(allhash),'new',new,'reuse',reuse)
