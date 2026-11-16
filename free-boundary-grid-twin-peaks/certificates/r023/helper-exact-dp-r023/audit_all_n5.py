import json,hashlib
from fractions import Fraction
from pathlib import Path
B=Path(__file__).parent;D=B/"n5-all-checkpoints"
src=hashlib.sha256((B/"helper_dp.cpp").read_bytes()).hexdigest()
fs=sorted(D.glob("*.json"));assert len(fs)==172
z=[json.loads(f.read_text()) for f in fs]
assert [x["record_index"] for x in z]==list(range(172))
assert all(x["source"]==src and x["cover_records"]==1 for x in z)
keys=[tuple(x["args"]) for x in z];assert len(set(keys))==172
N=25;num={}
for x in z:
 for k,v in x["rank_numerators"].items():num[k]=num.get(k,0)+int(v)
gain=sum((Fraction(v,N-int(k)-1) for k,v in num.items()),Fraction())
words=sum(int(x["cover_word_total"]) for x in z)
out={"status":"passed","records":172,"word_total":str(words),"gain_num":str(gain.numerator),"gain_den":str(gain.denominator),"avg_num":str((gain/words).numerator),"avg_den":str((gain/words).denominator),"source":src,"rank_numerators":{k:str(v) for k,v in num.items()}}
(B/"n5-all-audit.json").write_text(json.dumps(out,indent=2)+"\n");print(json.dumps(out,indent=2))
