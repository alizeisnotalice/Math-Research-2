#!/usr/bin/env python3
"""Audit the complete exact n=7 corner-monotone prefix scan receipt."""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path

HERE=Path(__file__).resolve().parent
R=HERE; RAW=HERE/'raw'
EXPECTED_SOURCE='850dea773006f6fccecaf3ded0980604071da056c01dd595ded5873a3f558f42'
EXPECTED_TOP=23386940526059745800104592446335589635840
EXPECTED_N6=20779014929894865426585600
def sha(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def resource(p:Path):
    s=p.read_text(); a={}
    first=s.splitlines()[0].split()
    a['real_seconds']=float(first[0])
    for ln in s.splitlines():
        if 'maximum resident set size' in ln:a['max_rss_bytes']=int(ln.split()[0])
        if ln.strip().endswith('swaps'):a['swaps']=int(ln.split()[0])
        if 'peak memory footprint' in ln:a['peak_memory_footprint_bytes']=int(ln.split()[0])
    return a

manifest=json.loads((R/'n7-monotone-manifest.json').read_text())
metadata=json.loads((R/'run-metadata.json').read_text())
assert manifest['n']==7 and manifest['path_count']==64 and len(manifest['rows'])==64
assert metadata['manifest_sha256']==sha(R/'n7-monotone-manifest.json')
assert metadata['source_sha256']==EXPECTED_SOURCE==sha(R/'completion_ie_selected.cpp')
bridge=json.loads((R/'smoke-n6-top.json').read_text())
assert int(bridge['hS'])==EXPECTED_N6
rows=[]; failures=[]
for m in manifest['rows']:
    w=m['word']; p=RAW/(w+'.json'); ep=RAW/(w+'.exit.json'); rp=RAW/(w+'.resource.txt'); sp=RAW/(w+'.started.json')
    d=json.loads(p.read_text()); e=json.loads(ep.read_text()); res=resource(rp)
    expected_cells=m['cells']; actual_cells=d['S_cells']
    checks={
      'word_matches':e['word']==w,
      'mask_matches':e['mask_hex']==m['mask_hex'],
      'cells_matches':actual_cells==expected_cells,
      'exit_zero':e['exit_code']==0,
      'output_hash_matches':sha(p)==e['output_sha256'],
      'resource_hash_matches':sha(rp)==e['resource_sha256'],
      'zero_swaps':res.get('swaps')==0,
      'n_is_7':d['n']==7,
    }
    if not all(checks.values()): failures.append({'word':w,'checks':checks})
    rows.append({'word':w,'mask_hex':m['mask_hex'],'hS':d['hS'],'independent_sets':d['independent_sets'],
                 'resource':res,'checks':checks,'output_sha256':sha(p),'resource_sha256':sha(rp)})
top=next(x for x in rows if x['word']=='RRRRRR')
values=[int(x['hS']) for x in rows]; min_h=min(values); max_h=max(values)
counter=[x['word'] for x in rows if int(x['hS'])<EXPECTED_TOP]
audit={'scope':'Exact finite h(S) scan for 64 corner-monotone n=7 paths; no asymptotic conclusion.',
       'comparison':'Only hS is compared. The program field candidate_holds is intentionally ignored.',
       'source_sha256':sha(R/'completion_ie_selected.cpp'),
       'manifest_sha256':sha(R/'n7-monotone-manifest.json'),'driver_exit':json.loads((R/'driver.exit.json').read_text()),
       'bridge_n6_hS':bridge['hS'],'top_hS':top['hS'],'expected_top_hS':str(EXPECTED_TOP),
       'min_hS':str(min_h),'minimizing_words':[x['word'] for x in rows if int(x['hS'])==min_h],
       'max_hS':str(max_h),'counterexample_words_hS_lt_top':counter,
       'all_hS_ge_top':not counter,'case_count':len(rows),'failed_case_checks':failures,
       'total_real_seconds':sum(x['resource']['real_seconds'] for x in rows),
       'max_rss_bytes':max(x['resource']['max_rss_bytes'] for x in rows),
       'rows':rows}
audit['all_checks_pass']=(not failures and audit['driver_exit']['exit_code']==0 and int(top['hS'])==EXPECTED_TOP and len(rows)==64)
(HERE/'audit-receipt.json').write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:audit[k] for k in ('all_checks_pass','case_count','top_hS','min_hS','minimizing_words','max_hS','counterexample_words_hS_lt_top','total_real_seconds','max_rss_bytes')},indent=2))
