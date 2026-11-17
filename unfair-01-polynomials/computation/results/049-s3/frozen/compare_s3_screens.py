"""Require complete independent replay of every necessary-screen input."""
import argparse,gzip,hashlib,json,math,time
from collections import Counter
from pathlib import Path

def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser()
    for n in ['parent','independent','parent-source','independent-source','prefix','out']:ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();start=time.monotonic();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();h=lambda b:hashlib.sha256(b).hexdigest()
    triples=[]
    for T in range(-16,17):
        for D in range(1,287):
            if (T*T-D)%4 or math.isqrt(D)**2==D:continue
            if 144*D<(203-12*abs(T))**2:
                N=(T*T-D)//4;assert -71<=N<=71;triples.append((T,N,D))
    triples.sort();assert len(triples)==669 and len(set(triples))==669
    counts=Counter();survivors=[];attest=[];stream=hashlib.sha256();prefix_hash=sha(a.prefix)
    for i,tn in enumerate(triples):
        p=a.parent/f'tn-{i:03}.json';q=a.independent/f'tn-{i:03}.json';x=json.loads(p.read_text());y=json.loads(q.read_text())
        assert x['status']==y['status']=='complete' and x['index']==y['index']==i and x['triple']==y['triple']==list(tn)
        assert x['binding']=={'prefix':prefix_hash,'source':sha(a.parent_source)}
        assert y['binding']=={'prefix':prefix_hash,'source':sha(a.independent_source)}
        for k in ['counts','survivors','raw_code_sha256']:assert x[k]==y[k],(i,k)
        packed=(a.parent/f'tn-{i:03}-codes.bin.gz').read_bytes();assert h(packed)==x['ledger_sha256']
        raw=gzip.decompress(packed);assert h(raw)==x['raw_code_sha256'] and len(raw)==x['counts']['f']
        freq=Counter(raw);assert set(freq)<=({0,255}|set(range(1,8))|set(range(21,35))|set(range(41,55)))
        assert freq[0]==x['counts'].get('prefix_rejected',0) and freq[255]==len(x['survivors'])==x['counts'].get('annulus_survived',0)
        for j in range(7):assert freq[1+j]==x['counts'].get(f'axis_{j}',0)
        assert sum(freq[z] for z in range(21,35))==x['counts'].get('outer_rejected',0)
        assert sum(freq[z] for z in range(41,55))==x['counts'].get('inner_rejected',0)
        for row in x['survivors']:assert row[:3]==list(tn)
        assert len({tuple(row) for row in x['survivors']})==len(x['survivors'])
        counts.update(x['counts']);survivors.extend(x['survivors']);stream.update(raw)
        attest.append({'index':i,'TN':tn,'parent_sha256':sha(p),'independent_sha256':sha(q),'raw_code_sha256':h(raw),'mothers':len(x['survivors'])})
    assert len({tuple(z) for z in survivors})==len(survivors)
    assert {k:counts[k] for k in ['TN','c','d','e','s7','f']}=={'TN':669,'c':8697,'d':213411,'e':3150761,'s7':2083109,'f':27909836}
    output={'status':'passed_all_669_complete_independent_screens_and_failure_bytes','counts':dict(counts),'surviving_mothers':len(survivors),'ordered_survivors_sha256':h((json.dumps(survivors,separators=(',',':'))+'\n').encode()),'all_failure_bytes_sha256':stream.hexdigest(),'per_TN':attest,'elapsed_seconds':time.monotonic()-start,'verifier_sha256':sha(Path(__file__))}
    a.out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({k:output[k] for k in ['status','counts','surviving_mothers','elapsed_seconds']}))

if __name__=='__main__':main()
