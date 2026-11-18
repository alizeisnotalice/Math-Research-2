"""Replay the scored original-leaf screen with signed rigid translations and angles."""
import argparse
import hashlib
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def compile_screen(base, original, main, work, cache):
    temp=work/'signed-source'
    temp.mkdir(exist_ok=True)
    v=(main/'five_grouped_localization.cpp').read_text()
    assert hashlib.sha256(v.encode()).hexdigest()==base.SOURCE_SHA
    def change(a,b):
        nonlocal v
        assert v.count(a)==1,a
        v=v.replace(a,b)
    change('        pts.clear();hull.clear();conv.clear();int off=0,ends[5]{};', '        L negative_angle[5]{},positive_angle[5]{},cos_width[5]{};\n        pts.clear();hull.clear();conv.clear();int off=0,ends[5]{};')
    token='theta[b]=small(up(std::max(std::abs(1024*lo[k]-anchor),std::abs(1024*hi[k]-anchor)),1024*ratio));'
    change(token,token+'''
                    negative_angle[b]=small(up(std::max<L>(0,anchor-1024*lo[k]),1024*ratio));
                    positive_angle[b]=small(up(std::max<L>(0,1024*hi[k]-anchor),1024*ratio));
                    cos_width[b]=small(up(W(theta[b])*theta[b],2*Q));''')
    change('W twice=0,error=0,grot','W rotation_remainder=0,twice=0,error=0,grot')
    t='error+=cross+std::min(dot,up(dot*(theta[k]+theta[j]),2*Q));'
    change(t,t+'\n                rotation_remainder+=cross+std::min(dot,up(dot*(theta[k]+theta[j]),2*Q));')
    change('        W best=-1;', (HERE/'tomar_signed_translation_screen.inc').read_text()+'\n        W best=-1;')
    (temp/'five_grouped_localization.cpp').write_text(v)
    base.SOURCE_SHA=hashlib.sha256(v.encode()).hexdigest()
    return original(temp,work,cache)

def main():
    if not __debug__:
        raise SystemExit("Assertions required.")
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("computation","tree","table","seeds","seed_receipt","work","cache"):
        parser.add_argument(name,type=Path)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--parts")
    parser.add_argument("--until",type=int,choices=(1,2,3,4),default=4)
    args=parser.parse_args()
    assert 1<=args.workers<=12
    sys.path.insert(0,str(args.computation.resolve()))
    import replay_tomar_083731_scored as base
    original=base.kernel
    base.kernel=lambda main,work,cache:compile_screen(base,original,main,work,cache)
    base.run(args)

if __name__=="__main__":
    main()
