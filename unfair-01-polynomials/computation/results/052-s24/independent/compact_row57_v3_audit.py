#!/usr/bin/env python3
"""Compact, per-field audit for the independent S={2,4} row-5/7 sieve.

It never builds a global state list.  Each field is checked against the
parent's independently implemented compact state encoding, then released.
"""
import hashlib,json,sys,time
from collections import Counter
from pathlib import Path

OWN=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/algebraic_global/round052-s24/row57-full-v3')
PARENT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round052-parent/row7-full')
OUT=OWN/'compact-audit.json'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def canon_own(r):
 return (r['a_m'],r['a_n'],r['h_m'],r['h_n'],r['g2'],r['q2'],r['g4'],r['q4'],
         r['R6'],r['S6'],r['c'],r['e'],r['eta5'],r['theta5'],r['g7'],r['q7'],
         tuple(r['f_values_passing_s6']))
def main():
 started=time.monotonic(); ownfiles=sorted(OWN.glob('field-*-D*.json'))
 assert len(ownfiles)==53, len(ownfiles)
 total=Counter(); reports=[]; totalrows=0
 for op in ownfiles:
  x=json.loads(op.read_text()); assert x['status']=='complete'
  D=x['binding']['D']; pp=PARENT/f'field-D{D}.json'; y=json.loads(pp.read_text())
  assert y['status']=='complete' and y['D']==D
  ox=[canon_own(r) for r in x['rows']]; py=[tuple(r[:-1])+ (tuple(r[-1]),) for r in y['rows']]
  os,ps=set(ox),set(py)
  same=(len(ox)==len(os)==len(py)==len(ps) and os==ps)
  reports.append({'D':D,'own_file':op.name,'own_sha256':sha(op),'parent_file':pp.name,
                  'parent_sha256':sha(pp),'own_rows':len(ox),'parent_rows':len(py),
                  'own_unique':len(os),'parent_unique':len(ps),'same_set_and_multiplicity':same,
                  'own_only_count':len(os-ps),'parent_only_count':len(ps-os)})
  if not same: raise AssertionError(reports[-1])
  total.update(x['counts']); totalrows+=len(ox)
  print(json.dumps({'D':D,'rows':len(ox),'same':same},sort_keys=True),flush=True)
 obj={'status':'complete','scope':'Per-field exact-state equality audit: independent row57 v3 versus parent compact row7 output. No global state list is materialized.',
      'own_source_sha256':sha(OWN.parent/'s24_row57_independent.py'),'parent_summary_sha256':sha(PARENT/'summary.json'),
      'fields':len(reports),'rows':totalrows,'counts':dict(total),'reports':reports,
      'python':sys.version,'elapsed_seconds':time.monotonic()-started}
 OUT.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':'complete','fields':len(reports),'rows':totalrows,'seconds':obj['elapsed_seconds']},sort_keys=True))
if __name__=='__main__':main()
