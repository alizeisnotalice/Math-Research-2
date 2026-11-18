"""Replay the common-fan AH chord with a free multiple at designated unconditioned roots."""
import argparse
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys

from certify_ah_deep_fan_chord import run as original, geometry, validate as inherited
from certify_ah_unconditional_root_fan import validate as free


def validate(ns,j,left,right,turn,sign,orientation,angular=12,depth=4):
    if orientation==3:
        return free(ns,j,left,right,turn,sign,1,angular=angular,depth=depth)
    if orientation==-3:
        return free(ns,j,left,right,turn,sign,-1,angular=angular,depth=depth)
    return inherited(ns,j,left,right,turn,sign,orientation,angular=angular,depth=depth)


def run(archive,inp,out):
    src=inspect.getsource(original)
    env={**original.__globals__,'validate':validate,'__file__':__file__}
    exec(compile(src,'breakpoint_fan_driver','exec'),env)
    env['run'](archive,inp,out)
    data=json.loads(out.read_bytes())
    for n in ('certify_ah_deep_fan_chord.py','certify_ah_unconditional_root_fan.py','certify_ah_affine_root_class.py'):
        data['source_hashes'][n]=sha256((archive/n).read_bytes()).hexdigest()
    out.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ('archive','input','output'):
        ap.add_argument(key,type=Path)
    args=ap.parse_args()
    run(args.archive.resolve(),args.input,args.output)
