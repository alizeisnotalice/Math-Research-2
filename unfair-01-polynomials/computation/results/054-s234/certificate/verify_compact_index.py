#!/usr/bin/env python3
"""Verify the compact S234 package without loading omitted full traces."""
from __future__ import annotations
import argparse, gzip, hashlib, json, os
from pathlib import Path

def sha(data): return hashlib.sha256(data).hexdigest()
def load_jsonl_gz(path, meta):
    blob=path.read_bytes(); assert len(blob)==meta['bytes'] and sha(blob)==meta['sha256']
    raw=gzip.decompress(blob); assert sha(raw)==meta['raw_sha256']
    return [json.loads(line) for line in raw.splitlines()]
def atomic(path,obj):
    temp=path.with_suffix(path.suffix+'.partial')
    with temp.open('w') as f: json.dump(obj,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(temp,path)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();root=a.package
    m=json.loads((root/'package-manifest.json').read_text())
    for x in m['copied_files']:
        p=root/x['path']; assert p.stat().st_size==x['bytes'] and sha(p.read_bytes())==x['sha256']
    fields=load_jsonl_gz(root/m['field_index']['path'],m['field_index'])
    records=[]
    for part in m['tree_index_parts']: records.extend(load_jsonl_gz(root/part['path'],part))
    assert len(fields)==m['counts']['fields']==58 and len(records)==m['counts']['root_survivors']==68593
    byD={}
    for r in records:
        assert r['status']=='dead_all_prefixes' and r['layers']==r['death_layer']>0 and r['max_live']>=0
        assert len(r['trace_sha256'])==64 and len(r['final_state_sha256'])==64
        byD.setdefault(r['D'],[]).append(tuple(r['row']))
    assert set(byD) <= {f['D'] for f in fields} and len({(r['D'],tuple(r['row'])) for r in records})==len(records)
    layers=sum(r['layers'] for r in records); maxdepth=max(r['death_layer'] for r in records); maxlive=max(r['max_live'] for r in records)
    assert (layers,maxdepth,maxlive)==(m['counts']['tree_layers'],m['counts']['max_depth'],m['counts']['max_live'])
    for f in fields:
        D=f['D']; ann=json.loads((root/f'input/annulus/field-D{D}.json').read_text())
        ab=json.loads((root/f'audit/ab/field-D{D}.json').read_text())
        assert [tuple(x) for x in ann['survivors']]==byD.get(D,[])
        assert ann['counts'].get('survived',0)==f['survivors']==len(byD.get(D,[]))==ab['mothers']
        assert ab['status']=='passed' and ab['counts'].get('dead_all_prefixes',0)==ab['mothers']
        assert ab['binding']==f['ab_binding'] and f['ab_binding']['HGQ_receipt']==f['hgq_field_sha256']
    out={'status':'passed_compact_s234_index','fields':len(fields),'records':len(records),'layers':layers,
         'max_depth':maxdepth,'max_live':maxlive,'manifest_sha256':sha((root/'package-manifest.json').read_bytes()),
         'scope':'Checks copied annulus/AB receipts and compact index coverage. Full HGQ traces remain intentionally external and are represented by per-record trace SHA-256.'}
    a.out.parent.mkdir(parents=True,exist_ok=True);atomic(a.out,out);print(json.dumps(out))
if __name__=='__main__':main()
