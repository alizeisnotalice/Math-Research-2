"""Compare all fresh 0.83750 packets to archived records after independent geometry replay."""
import hashlib,json,sys,zipfile
from pathlib import Path
P=Path(__file__).resolve().parent;fresh=Path(sys.argv[1]);sha=lambda b:hashlib.sha256(b).hexdigest()
old=(P/'records/stage-ledger.json').read_bytes();assert sha(old)=='3e83c009e5776e2af17be09d3d712d0cd5bc83236df439ed02cb83559fec1dd9';old=json.loads(old);new=json.loads((fresh/'stage-ledger.json').read_bytes())
assert {k:v for k,v in old.items() if k!='binary_sha'}=={k:v for k,v in new.items() if k!='binary_sha'}
pins=['b78fa27548a7fc95c269d1901418aa6016a5a2e2619e9ee3e1d073eeabcceb9e','d987d17b50ec3f37c04e9b6682952658096501943e60375cb3abb602ce066edc','5f87230a15f255984260750a43753bbcfe0ac0308ce616d05c4016bb95530437','902f6b22beac1a28a5e07236535285c7e26d570845db033dd8fb544bfbedfc0f']
for i,pin in enumerate(pins,1):
 path=P/f'volumes/screen{i}.zip';assert sha(path.read_bytes())==pin
 with zipfile.ZipFile(path) as z:
  rows=json.loads(z.read('MANIFEST.json'));names={f'stage{i}/part{part:04}.{ext}' for part,_ in old['audits'][i-1]['work_files'] for ext in ('npz','json')};assert len(rows)==len(names)==2048 and {r[0] for r in rows}==names and set(z.namelist())==names|{'MANIFEST.json'}
  for name,size,digest in rows:
   data=z.read(name);assert len(data)==size and sha(data)==digest and (fresh/name).read_bytes()==data, name
 a=(P/'records/retained.json').read_bytes();assert sha(a)=='06c7b4202f85d404c5e04cd53d19f4e4bd51a81da8e355e85f1f05a2913dcbc0';assert (fresh/'retained.json').read_bytes()==a
 print('verified layer',i,'files',len(rows),flush=True)
print(json.dumps(dict(complete=True,parts=1024,layers=4,remaining=5442,retained_sha256=sha(a),scope='fresh kernel output matches archived per-shard output; inherited seed geometry remains a separate input')))
