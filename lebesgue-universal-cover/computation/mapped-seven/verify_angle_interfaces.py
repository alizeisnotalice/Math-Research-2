"""Check signed rotations, the circular seam, and native mixed-grid sharing."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

from flint import arb,ctx
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_angle_rows import CircleMap,ball


def sha(data):return sha256(data).hexdigest()


def boundaries(src):
    p=(F(1,3),F(1,4))
    def r(q,a=F(4,5)):
        x,y=q
        return F(3,5)*x-a*y,a*x+F(3,5)*y
    def refl(q):return q[0],-q[1]
    h0=r(refl(p))
    incorrect=refl(r(p))
    assert h0!=incorrect and h0==refl(r(p,-F(4,5)))
    proof={}
    for file in ("P962.json","Q962.json"):
        raw=(src/file).read_bytes()
        m=json.loads(raw)
        b=[F(m["epsilon"],)/m["den"]*q for q in m["angles"]]
        maps=[CircleMap(b[i:i+384]) for i in range(0,len(b),384)]
        h={}
        for idx,v in enumerate(maps):
            assert v.shift(F(0))==v.shift(F(1))==v.shift(F(2))
            s=v.shift(F(1535,1536))
            assert v.shift(F(-1,1536))==s
            assert v.shift(F(1,1536))==v.shift(F(1537,1536))
            assert v.shift(F(383,384))==v.shift(F(4*383,1536))
            seam=v.shift(F(1))-v.shift(F(383,384))
            # Past the seam, the next coarse cell uses the same first knot.
            assert v.shift(F(1))+F(1,4)*(v.shift(F(385,384))-v.shift(F(1)))==v.shift(F(1537,1536))
            h[str(idx)]=dict(last_coarse_shift_delta=str(seam),first_fine_on_seam=str(s))
        if len(maps)==2:
            mismatches=[maps[1].shift(F(i,384))+maps[0].shift(-F(i,384)) for i in range(384)]
            assert any(mismatches)
            h["difference_if_wrongly_forced_reflected_map"]=str(max(map(abs,mismatches)))
        h["sha256"]=sha(raw)
        proof[file]=h
    slowed=F(".991")
    threshold=F(39,2)
    assert 2*arb.pi()**2*ball(slowed*slowed)<ball(threshold)<2*arb.pi()**2
    curvature_with_additive_shift=F(20)
    assert 2*arb.pi()**2<ball(curvature_with_additive_shift)<2*arb.pi()**2*ball("1.02")
    return dict(exact_rotation_then_reflection=[str(z) for z in incorrect],
                exact_reflection_then_rotation=[str(z) for z in h0],
                rejected_unsquared_speed_rule=dict(test_slope=str(slowed),desired_curvature=str(threshold)),
                rejected_virtual_in_curvature_rule=dict(physical_min="1",artificial_shift="1/50",
                                                        desired_curvature=str(curvature_with_additive_shift)),
                seams=proof)


def check_interfaces(comp, src):
    archives=json.loads((src/"source-interface.json").read_bytes())["archives"]
    kernels={}
    def inspect(z,prefix):
        for n in z.namelist():
            if n.endswith(".zip"):
                with ZipFile(BytesIO(z.read(n))) as child:inspect(child,prefix+"::"+n)
            elif n.endswith((".cpp",".hpp")):
                raw=z.read(n)
                if b'constexpr int CIRCLE=384' not in raw:
                    continue
                expected=(
                    b'const int offset[4]={0,CIRCLE,6*CIRCLE,16*CIRCLE};',
                    b'const int span[4]={CIRCLE,CIRCLE,2*CIRCLE,2*CIRCLE};',
                    b'const int refined[4]={-1,2*CIRCLE,8*CIRCLE,-1};',
                    b'if(!full_cell(refined[f]+4*angle+h))',
                )
                # R4 chooses only (whole physical 384) coarse groups and lacks a refined path.
                direct=b'const int span[4]={CIRCLE,CIRCLE,2*CIRCLE,2*CIRCLE};' not in raw
                if not direct:
                    assert all(b in raw for b in expected),n
                elif b'1536' in raw or b'refined[f]' in raw:
                    print('unclassified layout',prefix,n,flush=True)
                    raise AssertionError(n)
                if b'const L bow=' in raw:
                    assert b'const L bow=18*std::min(left,right)/(knots*knots);' in raw,n
                    assert b'const L knots=(i>=768 && i<2304) || (i>=3072 && i<6144) ? 1536:384;' in raw,n
                entry=kernels.setdefault(sha(raw),dict(files=[],use_fine=not direct,bow=b'const L bow=' in raw,
                    snippets=[x.strip() for x in raw.decode().splitlines() if any(p in x for p in ('int offset','int span','int refined','int CIRCLE','knots=','const L bow=','cells!=','4*angle','for(int f=0;'))]))
                entry["files"].append(prefix+"::"+n)
    for name,pin in archives.items():
        raw=(comp/name).read_bytes()
        assert sha(raw)==pin
        with ZipFile(BytesIO(raw)) as z:inspect(z,name)
    assert len(kernels)>=4 and any(r["bow"] for r in kernels.values())
    return kernels


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("comp","source","out"):
        p.add_argument(k,type=Path)
    a=p.parse_args()
    ctx.prec=256
    result=boundaries(a.source)
    kernels=check_interfaces(a.comp,a.source)
    result["interfaces"]=kernels
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("archived native kernels",len(kernels), "with shared bow",sum(x["bow"] for x in kernels.values()))
    print(json.dumps(result["seams"],indent=2))
