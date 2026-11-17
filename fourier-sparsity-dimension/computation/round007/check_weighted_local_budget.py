#!/usr/bin/env python3
"""Exact rational local-support budget on all stored n<=4 maximal code witnesses."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'round006'))
from verify_switch_closure import echelon,fwht


def span(gs):
    out=[0]
    for g in gs:out += [x^g for x in out]
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    raw=a.input.read_bytes();data=json.loads(raw);counts=Counter();examples=[]
    for per in data['by_n'].values():
        for row in per['rows']:
            counts['groups']+=1;S=row['support'];c=row['integercoeffs'];N=1<<row['n'];s=len(S)
            G=row['C_basis']+row['lifted_signcode_basis'];V=span(G);assert len(set(V))==1<<row['dimV']
            cache={}
            def partial(mask):
                if mask not in cache:
                    vec=[0]*N
                    for i,alpha in enumerate(S):
                        if mask>>i&1:vec[alpha]=c[i]
                    cache[mask]=fwht(vec)
                return cache[mask]
            flat=len(set(abs(x) for x in c))==1
            for T in V:
                if not T:continue
                counts['local_subspaces']+=1
                proper=T!=(1<<s)-1
                counts['proper_supports']+=proper
                counts['nonflat_proper_supports']+=proper and not flat
                local=[v for v in V if not(v&~T)]
                B=[x for x,_ in echelon(local).values()]
                assert len(local)==1<<len(B) and T in local
                atoms={}
                for i in range(s):
                    if T>>i&1:
                        label=sum(((b>>i)&1)<<j for j,b in enumerate(B))
                        assert label
                        atoms[label]=atoms.get(label,0)|(1<<i)
                energy=Fraction(sum(x*x for i,x in enumerate(c) if T>>i&1),N*N)
                bound=T.bit_count()*energy
                assert len(B)<=len(atoms)<=bound
                PT=partial(T);assert set(PT)<={-N,0,N}
                assert Fraction(sum(x!=0 for x in PT),N)==energy
                for atom in atoms.values():
                    counts['atoms']+=1
                    assert atom.bit_count()*energy>=1
                    counts['atom_uncertainty_equalities']+=atom.bit_count()*energy==1
                    PA=partial(atom)
                    assert all(not PA[x] for x in range(N) if not PT[x])
                if proper and not flat and len(B)>1 and len(examples)<5:
                    examples.append({'n':row['n'],'s':s,'support':S,'integercoeffs':c,'T_hex':hex(T),'dimG':len(B),'atoms':len(atoms),'energy':str(energy),'bound':str(bound)})
    result={'status':'passed','arithmetic':'integers and Fraction','tolerance':0,'counts':dict(counts),'examples':examples,'python':sys.version,'platform':platform.platform(),'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','counts':dict(counts)}))

if __name__=='__main__':main()
