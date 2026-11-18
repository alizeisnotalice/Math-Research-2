"""Certify open neighborhoods on the 948 diagonal by log-ratio interpolation."""
import argparse,json,os,time
from datetime import datetime,UTC
from pathlib import Path
from certify import D,LO,HI
from certify_midpoint import certificate


def verify(idx,out,gate):
    j,k,y,*_=json.loads(Path('diagonal-948.json').read_text())['weakest_refined'][idx]
    start_center=D(str(y)).quantize(D('.001'))
    if time.time()<gate:
        time.sleep(gate-time.time())
    started=time.time(); cpu=time.process_time()
    samples=[]
    for i in range(-6,7):
        x=start_center+D(i)*D('.001')
        data=certificate(j,k,x,x,12)
        p=D(data['retained_positive_scaled'][0]);n=D(data['negative_scaled'][1]);assert p>0
        b=LO.subtract(1,HI.divide(n,p))
        samples.append([str(x),str(b),data['retained_positive_scaled'],data['negative_scaled']])
    boxes=[]
    for left,right in zip(samples,samples[1:]):
        bound=LO.subtract(min(D(left[1]),D(right[1])),D('0.0000255'))
        boxes.append([left[0],right[0],str(bound),bound>0])
    out.write_text(json.dumps(dict(j=j,k=k,pid=os.getpid(),started=started,ended=time.time(),cpu=time.process_time()-cpu,
        points=samples,boxes=boxes,all_pass=all(b[3] for b in boxes)),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('shard',type=int);p.add_argument('out',type=Path);p.add_argument('--gate',type=float,default=datetime(2026,11,18,7,20,tzinfo=UTC).timestamp());a=p.parse_args();verify(a.shard,a.out,a.gate)
