from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib,json,subprocess,time

def sha(file):return hashlib.sha256(file.read_bytes()).hexdigest()

def main():
    p=Path(__file__).parent
    exe=p/'s3_integer_enclosure'
    src=p/'s3_integer_enclosure.cpp'
    assert sha(src)=='7710391531c9bfff04b807e70c95c87695ea0cbe90d15e44281bdb60c559d8df'
    index=json.loads((p/'manifest.json').read_text())
    assert index['groups']==1915327 and len(index['parts'])==13
    start=time.monotonic()
    def check(i):
        f=p/i['file'];assert f.stat().st_size==i['bytes'] and sha(f)==i['sha256']
        output=f.with_suffix('.pass.jsonl');log=f.with_suffix('.counts.json')
        with output.open('xb') as a,log.open('xb') as b:
            subprocess.run([str(exe),str(f)],stdin=subprocess.DEVNULL,stdout=a,stderr=b,check=True)
        counts=json.loads(log.read_bytes())
        with output.open('rb') as x:assert counts['unresolved']+counts['pass']==sum(1 for _ in x)
        return {'input':i,'counts':counts,'output':output.name,'output_sha256':sha(output)}
    with ThreadPoolExecutor(max_workers=13) as pool:
        rows=list(pool.map(check,index['parts']))
    keys=rows[0]['counts']
    counts={k:sum(r['counts'][k] for r in rows) for k in keys}
    assert counts['groups']==1915327 and counts['mothers']==1523605945
    summary={'status':'complete', 'workers':13, 'seconds':time.monotonic()-start,
             'manifest_sha256':sha(p/'manifest.json'),'source_sha256':sha(src),'binary_sha256':sha(exe),
             'counts':counts,'parts':rows,'child_processes_waited':True}
    (p/'summary.partial').write_text(json.dumps(summary,indent=2)+'\n');(p/'summary.partial').replace(p/'summary.json')
    print(json.dumps(counts),flush=True)
if __name__=='__main__': main()
