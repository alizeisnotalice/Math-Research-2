"""Exact support values of the compulsory edge ends on both near-rotated normals."""
from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import json
from flint import arb,ctx
def main():
    ctx.prec=192;deg=arb.pi()/180
    s=arb('1.52')*deg;d=arb('.02')*deg
    t=((arb.pi()/6-s)/2).tan()/2
    common=t*d.sin();loss=(1-d.cos())/2
    eplus=common-loss;eminus=common+loss
    assert eminus>eplus>0
    smax=arb('1.54')*deg
    tight=((arb.pi()/6-smax)/2).tan()/2
    assert tight>arb('.1267')
    # The source's lower bounds for the clockwise mandatory arcs hold on the
    # entire parent parameter interval [1.4,1.8] degrees.
    arcs={'TI':arb('.0118426'),'GK':arb('.0118750')}
    assert min(arcs.values())>arb('.04')*deg
    return dict(interval=['1.50','1.54'],compulsory_arc_reference='proofs/sigma-mandatory-arcs.md',
        coefficient_at_152=str(t),all_coefficient_lower=str(tight.lower()),
        sharp_offsets_150_154_at_152={'m3_clockwise':str(eminus),'m5_counterclockwise':str(eplus)},
        verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
if __name__=='__main__':
    p=ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args();r=main()
    a.out.write_text(json.dumps(r,indent=1)+'\n');print(r['sharp_offsets_150_154_at_152'])
