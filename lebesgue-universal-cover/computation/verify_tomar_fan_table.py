"""Freshly certify every stored integer witness, without the table-generation formula."""
from fractions import Fraction as F
import hashlib
from pathlib import Path
import numpy as np
from flint import arb,ctx

PINNED_MANIFEST="56cbe5dfdadc6b5aa07cf6ec9acf4f637a3c4fa78db54df0ae2d78916d292410"
def verified_table(path,certdir,header):
    ctx.prec=128
    q,grid,fine=2**28,1024,2**44
    manifest=(certdir/"SHA256SUMS").read_bytes()
    assert hashlib.sha256(manifest).hexdigest()==PINNED_MANIFEST
    with np.load(path,allow_pickle=False) as f:
        assert int(f["scale"])==q and int(f["grid"])==grid and int(f["original_scale"])==fine
        assert f["angle_indices"].tolist()==[2,5,8] and str(f["manifest_sha256"])==PINNED_MANIFEST
        for key,name in [("original_lo","lo"),("original_hi","hi")]:
            exact=[F(x)*fine for x in header[name]]
            assert all(x.denominator==1 for x in exact)
            assert f[key].dtype==np.int64 and f[key].tolist()==list(map(int,exact))
        points,counts=f["points"],f["counts"]
        generator_sha=str(f["generator_sha256"])
    assert generator_sha==hashlib.sha256(Path(__file__).with_name("tomar_fan_table.py").read_bytes()).hexdigest()
    assert counts.dtype==np.int32 and counts.shape==(5,) and (counts>=3).all() and (counts<=1024).all() and counts.sum()<=4096
    assert points.shape==(grid+1,int(counts.sum()),2) and points.dtype==np.int32
    assert (np.abs(points.astype(np.int64))<q).all()
    offsets=np.r_[0,np.cumsum(counts)]
    for body,n in enumerate((0,3,5,7,9)):
        block=points[:,offsets[body]:offsets[body+1],:]
        if body<2:
            assert (block==block[0]).all()
        angles=[0] if body<2 else range(grid+1)
        hi=F(0) if body<2 else F(header["hi"][3*body-4])
        for j in angles:
            theta=hi*j/grid;offset=arb(theta.numerator)/theta.denominator
            corners=[]
            if body:
                radius=1/(2*(arb.pi()/(2*n)).cos())
                corners=[(radius*(offset+2*arb.pi()*k/n).cos(),radius*(offset+2*arb.pi()*k/n).sin()) for k in range(n)]
            for uu,vv in block[j]:
                x,y=arb(int(uu))/q,arb(int(vv))/q
                if body:
                    assert all((x-X)**2+(y-Y)**2<1 for X,Y in corners)
                else:
                    assert x*x+y*y<arb(1)/4
    return points,counts
