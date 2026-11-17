"""Recompute all 669*13 prefix representatives through row six via HGQ."""
import argparse,hashlib,importlib.util,json,time
from pathlib import Path

def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--prefix',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();start=time.monotonic()
    spec=importlib.util.spec_from_file_location('hgq',a.source);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
    x=json.loads(a.prefix.read_text());triples=sorted({tuple(r[:3]) for r in x['per_TN_c']});assert len(triples)==669
    rows=[];counts={4:0,5:0,6:0}
    for T,N,D in triples:
        for c in range(-6,7):
            z=v.tree([T,N,D,c,0,0,0],depth=6)
            assert z['status'] in ('unknown_depth','dead_all_prefixes')
            for layer in z['trace']:
                if layer['k'] in counts:counts[layer['k']]+=layer['live']
            for s in z.get('live_prefixes',[]):
                G,Q=s['G'],s['Q'];rows.append([T,N,D,c,G[2]+c,Q[2]+c,*G,*Q])
    assert len(rows)==17760 and sorted(rows)==sorted(x['states'])
    assert counts=={4:26091,5:20070,6:17760}
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();out={'status':'passed_all_8697_representatives_and_17760_complete_prefix_states','counts':counts,'source_sha256':sha(a.source),'prefix_sha256':sha(a.prefix),'script_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
