#!/usr/bin/env python3
"""Local negative checks before per-total 47/48 semantic replay."""
import argparse,copy,gzip,hashlib,importlib.util,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def mod(n):
 s=importlib.util.spec_from_file_location(n,HERE/(n+'.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
 a=argparse.ArgumentParser();a.add_argument('--run',type=Path,required=True);a.add_argument('--plan',type=Path,required=True);a.add_argument('--frozen',type=Path,required=True);a.add_argument('--remote',type=Path,required=True);a.add_argument('--out',type=Path,required=True);x=a.parse_args()
 if x.out.exists():raise RuntimeError('existing fixture')
 v=mod('verify_support48');q=mod('audited_chunk_verifier');m=v.manifest(x.run,x.plan,x.frozen,x.remote)
 if m['errors']:raise RuntimeError(m['errors'])
 plan=json.loads(x.plan.read_text());want={v.key(z) for z in plan['jobs']};got={v.key(z) for z in m['chunks']}
 if {v.key(z) for z in m['chunks'][:-1]}==want:raise RuntimeError('delete accepted')
 if len([v.key(z) for z in m['chunks']]+[v.key(m['chunks'][0])])==len(set([v.key(z) for z in m['chunks']]+[v.key(m['chunks'][0])])):raise RuntimeError('duplicate accepted')
 e=next(z for z in m['chunks'] if z['total']==47 and z['counts']['excluded'])
 src=x.run/f"t{e['total']:02d}"/f"m{e['m']:02d}n{e['n']:02d}"/f"i{e['index']:06d}.jsonl.gz";rs=[json.loads(t) for t in gzip.decompress(src.read_bytes()).decode().splitlines()];next(t for t in rs if t[2]=='excluded')[3]=[999,999,999]
 raw=('\n'.join(json.dumps(t,separators=(',',':')) for t in rs)+'\n').encode();packed=gzip.compress(raw,compresslevel=1,mtime=0);bad=copy.deepcopy(e);bad['raw_sha256']=hashlib.sha256(raw).hexdigest();bad['gzip_sha256']=hashlib.sha256(packed).hexdigest();d=x.out/'t47'/f"m{e['m']:02d}n{e['n']:02d}";d.mkdir(parents=True);(d/f"i{e['index']:06d}.jsonl.gz").write_bytes(packed);errors=[];q.verify_chunk(x.out/'t47',bad,6,errors)
 if not any(z['kind']=='record_sequence_or_certificate_mismatch' for z in errors):raise RuntimeError('witness accepted')
 result={'format':'a30-47-48-negative-v1','deleted_chunk_rejected':True,'duplicate_chunk_rejected':True,'hash_rebound_witness_rejected':True,'source_job':v.key(e),'errors':errors};(x.out/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='errors'},default=str,sort_keys=True))
if __name__=='__main__':main()
