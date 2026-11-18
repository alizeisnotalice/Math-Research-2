"""Rank sampled five-body material chains by a floating-point first-variation area score."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import numpy as np

HERE=Path(__file__).resolve().parent
def sha(path):return sha256(path.read_bytes()).hexdigest()
def cpu_scores(center, width, moves):
    phi=(np.arange(64)+0.31783797)*2*np.pi/64
    u=np.stack((np.cos(phi),np.sin(phi)),axis=-1)
    values=[]
    for move in moves:
        point=.5*u.copy(); material=point.copy();owner=np.zeros(64,dtype=int);h=np.full(64,.5)
        for b,n in enumerate((3,5,7,9),1):
            j=3*(b-1); a0=0 if b==1 else center[j-1]
            da=0 if b==1 else width[j-1]*move[j-1]
            ds=2*np.pi/n;R=.5/np.cos(ds/4);angle=a0+da
            v=angle+ds*np.rint((phi-angle)/ds)
            direct=abs(phi-v)<=ds/4
            w=np.where(direct,v,angle+ds*np.rint((phi+np.pi-angle)/ds))
            a=R*np.stack((np.cos(w),np.sin(w)),axis=-1)+(~direct)[:,None]*u
            t=center[j:j+2]; shift=width[j:j+2]*move[j:j+2]
            val=np.sum((a+t+shift)*u,axis=1);mask=val>h
            ca,sa=np.cos(da),np.sin(da)
            ref=a@[[ca,-sa],[sa,ca]]
            point[mask]=t+ref[mask];material[mask]=ref[mask];owner[mask]=b;h[mask]=val[mask]
        nx=np.roll(point,-1,axis=0);dc=nx-np.roll(point,1,axis=0)
        value=np.sum(point[:,0]*nx[:,1]-point[:,1]*nx[:,0]);grad=np.zeros(11)
        for b in range(1,5):
            j=3*(b-1);mask=owner==b
            grad[j]=dc[mask,1].sum();grad[j+1]=-dc[mask,0].sum()
            if b>1:grad[j-1]=-np.sum(material[mask]*dc[mask])
        values.append(.5*(value-abs(grad)@width))
    return np.asarray(values)
def main():
    p=ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('metadata',type=Path);p.add_argument('work',type=Path)
    p.add_argument('--start',type=int,default=0);p.add_argument('--stop',type=int)
    p.add_argument('--native',type=Path,required=True)
    a=p.parse_args();meta=json.loads(a.metadata.read_bytes());assert sha(a.input)==meta['file_sha256']
    assert meta['coordinate_denominator']==2**44 and meta['original_box']
    stop=meta['rows'] if a.stop is None else a.stop
    assert 0<=a.start<stop<=meta['rows'];a.work.mkdir(exist_ok=True,parents=True)
    with np.load(a.input,allow_pickle=False) as z:
        lo=z['lo'][a.start:stop];hi=z['hi'][a.start:stop]
        part=z['part'][a.start:stop];pos=z['pos'][a.start:stop]
    center=(lo+hi)/2**45;width=(hi-lo)/2**45
    poses=np.hstack((center,width)).astype('<f4')
    rng=np.random.default_rng(20261118)
    corners=rng.integers(0,2,(96,11),dtype=np.int32)*2-1
    moves=np.vstack((np.zeros((1,11),np.int32),corners[:8])).astype('<i4')
    ids=np.array([0,*range(23,31)],np.uint8)
    x=a.work/'poses.f32';m=a.work/'moves.i32';output=a.work/'scores.f32'
    poses.tofile(x);moves.tofile(m)
    native=json.loads(subprocess.check_output([a.native.resolve(),HERE/'tomar_witness_rank.metal',x,m,output],text=True))
    assert native['rows']==len(part) and native['moves']==len(moves) and native['directions']==64
    values=np.fromfile(output,dtype='<f4').reshape(len(part),len(moves))
    assert np.all(np.isfinite(values))
    samples=np.unique(np.linspace(0,len(part)-1,min(11,len(part)),dtype=int))
    deltas=[]
    for i in samples:
        ref=cpu_scores(center[i],width[i],moves)
        deltas.append(float(np.max(np.abs(ref-values[i]))))
    # A near tie between two supporting bodies can change the first derivative
    # of a proposed material chain. The independently measured deviations are
    # recorded explicitly; the ranking is a numerical proposal.
    order=np.argsort(-values,axis=1,kind='stable')
    path=a.work/'witness-rank.npz'
    np.savez_compressed(path,part=part,pos=pos,move_order=ids[order],scores=values,move_ids=ids)
    report=dict(complete=True,status='floating witness proposal',start=a.start,stop=stop,native=native,
                input_sha256=sha(a.input),source_sha256=sha(HERE/'tomar_witness_rank.metal'),
                moves_sha256=sha(m),ranking_sha256=sha(path),float64_reference_rows=samples.tolist(),
                maximum_absolute_differences=deltas,reference_error_max=max(deltas),
                relative_to_center_improvement_count=int((np.max(values,axis=1)>values[:,0]+1e-7).sum()))
    (a.work/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,separators=(',',':')))
if __name__=='__main__':main()
