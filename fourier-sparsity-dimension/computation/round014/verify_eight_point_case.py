#!/usr/bin/env python3
"""Low-memory witness check and independent direct-sum upper certificate for three eight-point cases."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
def canonical(vs):
    piv={}
    for x in vs:
        for j in sorted(piv,reverse=True):
            if x>>j&1:x^=piv[j]
        if x:
            j=x.bit_length()-1
            for k in piv:
                if piv[k]>>j&1:piv[k]^=x
            piv[j]=x
    return list(piv.values())
def main():
    a=argparse.ArgumentParser()
    a.add_argument('--inputs',type=Path,required=True);a.add_argument('--index',type=int,required=True)
    a.add_argument('--result',type=Path,required=True);a.add_argument('--serial-binary',type=Path,required=True);a.add_argument('--parallel-binary',type=Path,required=True);a.add_argument('--out',type=Path,required=True)
    z=a.parse_args();task=json.loads(z.inputs.read_text())[z.index];result=json.loads(z.result.read_text())
    b=task['basis'];d=len(b);assert len(canonical(b))==d==result['D']
    union=0;terms=[]
    for text in result['partition']:
        cell=int(text,16);assert cell and not union&cell and len(canonical(b+[cell]))==d
        union|=cell
        dim=d-len(canonical([g&~cell for g in b]));terms.append(dim+1)
    assert union==(1<<64)-1 and sum(terms)==result['kappa']
    groups=[[[0,1,2,3],[4],[5],[6],[7]],[[0,1,2,3],[4,5],[6],[7]],[[0,1,2,3,4],[5],[6],[7]]][z.index]
    components=[]
    for A in groups:
        positions=[p for p in range(32) if p%8 in A]
        indicator=sum(1<<p for p in range(64) if p%8 in A)
        assert len(canonical(b+[indicator]))==d
        # Algebraic direct-sum test: every projected generator lies in the full code.
        assert all(len(canonical(b+[g&indicator]))==d for g in b)
        cb=canonical([sum(((g>>positions[i%len(positions)])&1)<<i for i in range(64)) for g in b])
        serial=subprocess.run([str(z.serial_binary.resolve()),str(len(cb)),*map(str,cb)],capture_output=True,text=True,check=True,timeout=10)
        parallel=subprocess.run([str(z.parallel_binary.resolve()),'4',str(len(cb)),*map(str,cb)],capture_output=True,text=True,check=True,timeout=10)
        s=json.loads(serial.stdout);p=json.loads(parallel.stdout);assert s['kappa']==p['kappa']
        components.append({'old_points':A,'D':len(cb),'kappa':s['kappa'],'serial_exit':serial.returncode,'parallel_exit':parallel.returncode})
    assert sum(x['kappa'] for x in components)==result['kappa']==24
    out={'status':'passed','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(z.inputs.read_bytes()).hexdigest(),'result_sha256':hashlib.sha256(z.result.read_bytes()).hexdigest(),'task':task['name'],'result':result,'partition_terms':terms,'direct_sum_components':components,'arithmetic':'F2 rank, integer DP','python':sys.version,'scope':'Witness lower bound and independent direct-sum upper certificate; not a general addback theorem.'}
    z.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','task':task['name'],'kappa':result['kappa'],'components':components}))
if __name__=='__main__':main()
