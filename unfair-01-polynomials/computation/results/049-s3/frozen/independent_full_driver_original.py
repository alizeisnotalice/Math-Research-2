import gzip,hashlib,json,os,time
from collections import Counter
from pathlib import Path
import independent_s3_audit as a
OUT=Path(__file__).parent/'full669';PARENT=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round048-parent/s3-screen-complete');LIMIT=14400
def atom(p,b):
 q=p.with_suffix(p.suffix+'.partial');q.write_bytes(b);os.replace(q,p)
def main():
 OUT.mkdir(exist_ok=True);pref=json.loads(a.PREFIX.read_text());allow={tuple(x[:4]) for x in pref['states']};all=a.ts();bind={'source':hashlib.sha256(Path(a.__file__).read_bytes()).hexdigest(),'prefix':hashlib.sha256(a.PREFIX.read_bytes()).hexdigest()};start=time.monotonic();done=[]
 for i,tr in enumerate(all):
  p=OUT/f'tn-{i:03}.json'
  if p.exists():
   x=json.loads(p.read_text());
   if x.get('binding')==bind and x.get('status')=='complete':done.append(i);continue
  if time.monotonic()-start>LIMIT:break
  c,b,s=a.one(tr,allow);x={'status':'complete','index':i,'triple':tr,'binding':bind,'counts':c,'survivors':s,'raw_code_sha256':hashlib.sha256(b).hexdigest(),'gzip_code_sha256':hashlib.sha256(gzip.compress(b,mtime=0)).hexdigest()}
  pp=PARENT/f'tn-{i:03}.json'
  if pp.exists():
   y=json.loads(pp.read_text())
   if y['counts']!=c or y['survivors']!=s or y['raw_code_sha256']!=x['raw_code_sha256']:
    atom(OUT/'MISMATCH.json',json.dumps({'index':i,'independent':x,'parent':y},sort_keys=True).encode());raise SystemExit(2)
  atom(p,(json.dumps(x,sort_keys=True)+'\n').encode());done.append(i)
 atom(OUT/'summary.json',(json.dumps({'status':'complete' if len(done)==669 else 'partial_unknown','completed':len(done),'indices':done,'binding':bind,'seconds':time.monotonic()-start},sort_keys=True)+'\n').encode())
if __name__=='__main__':main()
