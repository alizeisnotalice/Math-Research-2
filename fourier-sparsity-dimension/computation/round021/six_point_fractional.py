from pathlib import Path
import itertools,json,time,hashlib
start=time.monotonic();B=[3,7,11,12,13,14]
rows=[tuple(1-2*((x&b).bit_count()&1) for b in B) for x in range(16)]
allowed=[]
for a in itertools.product(range(-4,5),repeat=6):
 if all(sum(s*v for s,v in zip(row,a)) in (-16,0,16) for row in rows):allowed.append(a)
nontrivial=[a for a in allowed if any(a)];fractional=[a for a in nontrivial if any(abs(x)!=4 for x in a)]
out={'grid_denominator':4,'grid_points':9**6,'allowed_count':len(allowed),'nonzero_count':len(nontrivial),'nonsaturated_count':len(fractional),'allowed_numerators':allowed,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
Path(__file__).with_name('six-point-fractional.json').write_text(json.dumps(out,indent=2)+'\n');print({k:out[k] for k in ['allowed_count','nonzero_count','nonsaturated_count','elapsed_seconds']})
