#!/usr/bin/env python3
"""Portable CPU regression of the launched GPU52 inverse-transform formula."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--original-dir',type=Path,required=True);ap.add_argument('--continuation-source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    data=json.loads((a.original_dir/'inputs.json').read_text())
    original=a.original_dir/'search.wl'
    for source in [original,a.continuation_source]:
        text=source.read_text();assert 'ib=bvI.dc-bvR.ds' in text and 'ib=bvI.dc-avR.ds' not in text
    dc,ds=np.asarray(data['DC']),np.asarray(data['DS']);size=dc.shape[1];rows=[]
    for taskid,mask in [(0,525312),(1,1),(2,0),(2,23166975),(2,23166976),(3,17)]:
        task=data['tasks'][taskid];bits=np.array([0]+list(map(int,format(mask,f"0{task['groups']-1}b"))))
        ar=bits@np.asarray(task['LR']);ai=bits@np.asarray(task['LI']);br=(1-bits)@np.asarray(task['LR']);bi=(1-bits)@np.asarray(task['LI'])
        avr,avi=np.exp(ar)*np.cos(ai),np.exp(ar)*np.sin(ai);bvr,bvi=np.exp(br)*np.cos(bi),np.exp(br)*np.sin(bi)
        ac=avr@dc+avi@ds;bc=bvr@dc+bvi@ds;imA=avi@dc-avr@ds;imB=bvi@dc-bvr@ds
        real=[np.array([1.]),np.array([1.])]
        for bit,group in zip(bits,task['roots']):
            roots=np.array([complex(*z) for z in group]);factor=np.polynomial.polynomial.polyfromroots(roots);factor=factor/factor[0]
            assert np.max(np.abs(factor.imag))<1e-8
            j=0 if bit else 1;real[j]=np.polynomial.polynomial.polymul(real[j],factor.real)
        ea,eb=[np.pad(x,(0,size-len(x))) for x in real]
        error=float(max(np.max(abs(ac-ea)),np.max(abs(bc-eb)),np.max(abs(imA)),np.max(abs(imB))))
        assert error<3e-4
        rows.append(dict(task=taskid,mask=mask,max_reconstruction_or_imag_error=error))
    report=dict(status='passed_portable_cpu_formula_regression',numpy=np.__version__,samples=rows,max_error=max(x['max_reconstruction_or_imag_error'] for x in rows),hashes={"inputs":hashlib.sha256((a.original_dir/'inputs.json').read_bytes()).hexdigest(),"original_source":hashlib.sha256(original.read_bytes()).hexdigest(),"continuation_source":hashlib.sha256(a.continuation_source.read_bytes()).hexdigest(),"verifier":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},scope='Six numerical reconstruction samples plus exact source-expression checks; not complete mask coverage or an exact polynomial certificate.')
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
