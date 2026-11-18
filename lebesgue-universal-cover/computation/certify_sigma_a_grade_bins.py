"""Check real double hits in the A-grade candidate with an extended partner graph."""
import argparse,inspect,sys
from pathlib import Path
from certify_sigma_extended_graph import geometry
from certify_sigma_reachable_roots import run as ordinary


def run(archive,inp,out):
    env={**ordinary.__globals__,'geometry':geometry,'__file__':__file__}
    exec(compile(inspect.getsource(ordinary),'sigma_A_grade_real_roots','exec'),env);env['run'](archive,inp,out)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.input,opt.output)
