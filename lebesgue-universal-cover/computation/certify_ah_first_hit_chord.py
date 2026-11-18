"""Global first-hit normalized cap certificate with orientation-specific deficit budgets."""
import argparse,inspect,json,sys
from hashlib import sha256
from pathlib import Path
from certify_ah_chain_chord import run as original
from certify_ah_first_hit_chain import geometry,validate


def run(archive,inp,out):
    budgets={}
    def chart(*args,**kw):
        ns,code,aux=geometry(*args,**kw);budgets.update({j:[str(row[j]) for row in ns['labelled_deficits']] for j in range(len(ns['deficit_by_bin'])) if not ns['deficit_by_bin'][j].is_zero()});return ns,code,aux
    source=inspect.getsource(original);env={**original.__globals__,'geometry':chart,'validate':validate,'__file__':__file__};exec(compile(source,'first_hit_complete','exec'),env);env['run'](archive,inp,out)
    data=json.loads(out.read_bytes());data['first_hit_deficits_by_straight_or_paired_label']=budgets
    for name in ('certify_ah_first_hit_chain.py','certify_ah_centered_deficit.py'):data['source_hashes'][name]=sha256((archive/name).read_bytes()).hexdigest()
    out.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
