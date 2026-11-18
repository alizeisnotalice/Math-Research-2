"""Extract and verify the portable 0.83750 seeds to a new replay input directory."""
import hashlib,json,sys,zipfile
from pathlib import Path
P=Path(__file__).resolve().parent;out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False);sha=lambda b:hashlib.sha256(b).hexdigest();archive=P/'volumes/seed50.zip'
assert sha(archive.read_bytes())=='9738fff07eaed42db51bf27839e4ed76c457ade80d7a1e5352efb9f0e5c31905'
with zipfile.ZipFile(archive) as z:
 rows=json.loads(z.read('MANIFEST.json'))
 assert len(rows)==1026 and len({r[0] for r in rows})==1026 and set(z.namelist())=={r[0] for r in rows}|{'MANIFEST.json'}
 for name,size,digest in rows:
  p=Path(name);assert name in ('seed.json','full-index-report.json') or (len(p.parts)==2 and p.parts[0]=='seeds' and p.name.startswith('part') and p.suffix=='.npz');data=z.read(name);assert len(data)==size and sha(data)==digest;(out/p.parent).mkdir(exist_ok=True);(out/p).write_bytes(data)
assert sha((out/'seed.json').read_bytes())=='7c5790b0a38585d36d38b2d667943487fbb3767ce44de199950d1e4723bf3d00';s=json.loads((out/'seed.json').read_bytes());assert s['remaining']==699060 and len(s['rows'])==1024 and s['target']=='67/80'
print(json.dumps(dict(seed_parts=1024,remaining=s['remaining'],seed_sha256=sha((out/'seed.json').read_bytes()))))
