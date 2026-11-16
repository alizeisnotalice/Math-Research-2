"""Reproduce one necessary-support ledger; discovery is not its verifier."""
import argparse,gzip,json,math,os,time
from pathlib import Path
from propagation import propagate

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--p',type=int,required=True);ap.add_argument('--q',type=int,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data') and 1<a.p<a.q and a.q>=2*a.p and math.gcd(a.p,a.q)==1
    p,q=a.p,a.q;selected={0};total=count=unknown=0;start=time.monotonic()
    def walk(x):
        if x==q:yield tuple(sorted(selected));return
        if x-p not in selected and x-(q-p) not in selected:
            selected.add(x);yield from walk(x+1);selected.remove(x)
        yield from walk(x+1)
    with a.out.open('xb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as f:
        f.write(b'{"rows":[')
        for T in walk(1):
            total+=1;K=sorted(set(T)|{j+q for j in T});res=[0]*p
            for j in K:res[j%p]+=1
            if 1 in res:continue
            row=propagate(p,q,K)
            if count:f.write(b',')
            f.write(json.dumps(row,separators=(',',':')).encode());count+=1;unknown+=bool(row['unknown_indices'])
        metadata={'schema':'remote-q-p-propagation-v1','p':p,'q':q,'total_T':total,'count':count,'unknown':unknown,'all_unit':count-unknown,'elapsed_seconds':time.monotonic()-start}
        f.write(b'],'+json.dumps(metadata,separators=(',',':')).encode()[1:])
    print(json.dumps(metadata))
if __name__=='__main__':main()
