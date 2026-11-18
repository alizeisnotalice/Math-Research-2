"""Run integer hinge proofs on a user-specified scratch directory."""
import argparse,json,subprocess,sys
from hashlib import sha256
from fractions import Fraction as F
from pathlib import Path
from hinge_bank import bank,host,S,write_poly
from pose_boxes import make

def split(box):
    weights=[1000,1000,163,1000,1000,115,1000,1000,89,1000,1000]
    j=max(range(11),key=lambda i:weights[i]*(box[i][1]-box[i][0]))
    a,b=box[j];m=(a+b)//2
    return [box[:j]+[[a,m]]+box[j+1:],box[:j]+[[m,b]]+box[j+1:]]

def subboxes(path):
    data=json.loads(path.read_bytes());leaves=set(data['leaves']);assert len(leaves)==len(data['leaves'])
    out=[]
    def walk(addr,box):
        if addr in leaves:out.append(box);return
        assert any(v.startswith(addr) for v in leaves)
        for i,b in enumerate(split(box)):walk(addr+str(i),b)
    walk('',data['box']);assert len(out)==len(leaves)
    return out

def run(c,root,scratch,kernel,part):
    scratch.mkdir(parents=True,exist_ok=True);result={}
    raw=(root/'new-points.json').read_bytes()
    assert sha256(raw).hexdigest()=='45cb554f79d2c30f23cdc526f7cc1a333b6f98c3dbaa60234862d449bb66de05'
    pts=json.loads(raw)
    assert len(pts['points'])==11 and all((x-u)**2+(y-v)**2<=pts['denominator']**2 for x,y in pts['points'] for u,v in pts['points'])
    cases={
        'anchors':[(4096,837650,None),(4096,837750,[1,2])],
        'residual':[(4096,837500,None)],'tube':[(1024,837550,None),(4096,837725,None)],
        'subtree114819':[(4096,837500,None)], 'reverse':[(4096,838750,None)], 'star-fine':[(16384,837680,None)]}
    for name in part or cases:
        tree=name=='subtree114819'
        if tree:
            from corner_core import core
            rows=subboxes(root/'subtree114819.json')
            polygons=[core(c,box,64) for box in rows]
        elif name=='star-fine':
            rows=make(c,root,'anchors')[:1]
            polygons=[host(c,rows[0],384)]
        elif name=='reverse':
            pose=json.loads((root/'new-six-pose.json').read_bytes())
            rows=[[[v,v] for v in pose['center']]]
            polygons=[host(c,rows[0],192)]
        else:
            rows=make(c,root,name)
            polygons=[host(c,box,192) for box in rows]
        for count,target,indices in cases[name]:
            shape='p30' if name=='reverse' else 'new11'
            pb=scratch/f'angles-{shape}-{count}.dat'
            if (shape,count) not in run.banks:
                if shape=='p30':
                    src=c/'finite-point-leaf-input.json'
                    assert sha256(src.read_bytes()).hexdigest()=='b99b3bf66adf904b8b99ba1e1c7f8e75abc1ea6d8e00728b9553c20b961ffebc'
                    bank(src,count,pb,F('.45'),F('.32'),F('.55'))
                else:
                    bank(root/'new-points.json',count,pb)
                run.banks.add((shape,count))
            ids=range(len(rows)) if indices is None else indices
            ph=scratch/f'{name}-{target}.dat'
            with ph.open('w') as f:
                f.write(f'{S} {len(ids)}\n')
                for i in ids:write_poly(f,polygons[i])
            output=subprocess.check_output([kernel,ph,pb,str(target)],text=True)
            got=[list(map(int,line.split())) for line in output.splitlines()]
            assert got==[[i,2*count,2*count,-1] for i in range(len(ids))],(name,got)
            result[f'{name}-{target}']=dict(pose_boxes=len(ids),orientation_intervals=2*count,strict_target=str(F(target,1000000)),closed_cases=2*count*len(ids))
            print(name,result[f'{name}-{target}'],flush=True)
    return result
run.banks=set()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('computation',type=Path)
    p.add_argument('scratch',type=Path);p.add_argument('--kernel',type=Path,required=True)
    p.add_argument('--only',nargs='+',choices=['anchors','residual','tube','subtree114819','reverse','star-fine']);p.add_argument('--output',type=Path)
    a=p.parse_args();assert __debug__
    r=run(a.computation,Path(__file__).resolve().parent,a.scratch,a.kernel.resolve(),a.only)
    if a.output:a.output.write_text(json.dumps(r,indent=2)+'\n')
