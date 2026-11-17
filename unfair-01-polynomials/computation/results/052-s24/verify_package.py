"""Full S24 replay: complete lattice, prefix cover, root codes and cofactor trees."""
import argparse,concurrent.futures as cf,gzip,hashlib,importlib.util,json,os,shutil,sys,tempfile,time
from pathlib import Path
from collections import Counter,defaultdict
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(p,n):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def init(package):
    global ROW7,ROOTCHECK,TREE,AB,PACKAGE
    PACKAGE=Path(package);f=PACKAGE/'frozen';ROW7=module(f/'s24_row7.py','r7');ROOTCHECK=module(f/'verify_s24_annulus.py','rootcheck');TREE=module(f/'s24_hgq.py','tree');AB=module(f/'ab_prefix_bfs.py','ab')
def field(D,rows,records,tmp):
    tmp=Path(tmp);inputdir=tmp/'row7';rootdir=tmp/'root';auditdir=tmp/'root-audits';counts,states=ROW7.one(rows)
    source=inputdir/f'field-D{D}.json';source.write_text(json.dumps({'status':'complete','D':D,'rows':states},separators=(',',':'))+'\n');state_count=len(states);del states
    r=ROOTCHECK.one(D,str(inputdir),str(rootdir),str(auditdir),time.time()+14400,sha(PACKAGE/'frozen/verify_s24_annulus.py'));assert r['status']=='passed'
    doc=json.loads((rootdir/f'field-D{D}.json').read_text());assert [z['row'] for z in records]==doc['survivors']
    stats=Counter();layers=0;deadmax=0;maxlive=0
    for z in records:
        if z['status']=='dead_all_prefixes':k=z['death'];assert 1<=k<=92
        else:assert z['status']=='unknown_depth' and z['row']==[13,1,1,-2,1,3,-1,0] and z['death'] is None;k=96
        t=TREE.tree(z['row'],depth=k,seconds=60)
        assert t['status']==z['status'] and t['depth_reached']==k
        direct=AB.run_core(z['row'],depth=k,cap=25000,seconds=60)
        assert direct['status']==('dead' if z['status']=='dead_all_prefixes' else z['status'])
        assert direct['depth_reached']==k and direct['levels']==t['trace']
        digest=hashlib.sha256(json.dumps(t['trace'],separators=(',',':'),ensure_ascii=False).encode()).hexdigest();assert digest==z['trace_sha256']
        ml=max(v['live'] for v in t['trace']);assert ml==z['maxlive'];maxlive=max(maxlive,ml);layers+=k;stats[z['status']]+=1
        if z['status']=='dead_all_prefixes':deadmax=max(deadmax,k)
    source.unlink()
    return {'D':D,'row7_states':state_count,'prefix_counts':counts,'root_counts':r['counts'],'tree_stats':dict(stats),'layers':layers,'deadmax':deadmax,'maxlive':maxlive}
def main():
    if not __debug__:raise RuntimeError('no -O')
    ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert 1<=a.workers<=12;r=a.package;started=time.monotonic()
    if a.out.exists():raise FileExistsError('Refusing existing output; select a fresh path to avoid stale success receipts')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((r/'core-manifest.json').read_text())
    for z in manifest['files']:p=r/z['path'];assert p.stat().st_size==z['bytes'] and sha(p)==z['sha256']
    cert=r/'certificate';meta=json.loads((cert/'manifest.json').read_text());fields=meta['fields'];assert len(fields)==len(set(fields))==53
    records=defaultdict(list)
    for z in meta['record_parts']:
        p=cert/z['file'];assert sha(p)==z['sha256'];raw=gzip.decompress(p.read_bytes());assert hashlib.sha256(raw).hexdigest()==z['uncompressed_sha256'];n=0
        for line in raw.splitlines():v=json.loads(line);records[v['row'][0]].append(v);n+=1
        assert n==z['records']
    assert sum(map(len,records.values()))==29227 and len({tuple(v['row']) for seq in records.values() for v in seq})==29227
    with tempfile.TemporaryDirectory(prefix='s24-full-',dir=a.out.parent) as tmp:
        tmp=Path(tmp);f=r/'frozen';pair=module(f/'s24_low.py','pair');six=module(f/'s24_row6.py','six')
        low=pair.main(tmp/'low.json');assert low['counts']['raw_pairs']==774100 and len(low['rows'])==33575
        frozen=json.loads(gzip.decompress((r/'audits/s24-low.json.gz').read_bytes()));assert low['counts']==frozen['counts'] and low['rows']==frozen['rows']
        low6=six.main(tmp/'low.json',tmp/'six.json');assert len(low6['rows'])==20849
        frozen=json.loads(gzip.decompress((r/'audits/s24-row6.json.gz').read_bytes()));assert low6['counts']==frozen['counts'] and low6['rows']==frozen['rows']
        groups=defaultdict(list)
        for z in low6['rows']:groups[z['D']].append(z)
        assert sorted(groups)==fields
        for name in ['row7','root','root-audits']:(tmp/name).mkdir()
        for z in meta['fields_manifest']:
            p=cert/'frozen-rootdocs'/z['rootdoc'];assert sha(p)==z['rootdoc_sha256'];shutil.copy2(p,tmp/'root'/p.name)
            p=cert/'frozen-first-failure-codes'/z['codes'];assert sha(p)==z['codes_gzip_sha256'];shutil.copy2(p,tmp/'root'/p.name)
        prefix_counts=Counter();root_counts=Counter();tree_counts=Counter();layers=0;maxlive=0;maxdead=0;done=[]
        with cf.ProcessPoolExecutor(max_workers=a.workers,initializer=init,initargs=(str(r),)) as pool:
            futures=[pool.submit(field,D,groups[D],records[D],str(tmp)) for D in fields]
            for fut in cf.as_completed(futures):
                x=fut.result();done.append(x['D']);prefix_counts.update(x['prefix_counts']);root_counts.update(x['root_counts']);tree_counts.update(x['tree_stats']);layers+=x['layers'];maxlive=max(maxlive,x['maxlive']);maxdead=max(maxdead,x['deadmax']);print(json.dumps({'verified_fields':len(done),'D':x['D'],'tree_counts':dict(tree_counts)}),flush=True)
        assert sorted(done)==fields and prefix_counts['row7_states']==1796814 and root_counts['unique_mothers']==20209730 and root_counts['survived']==29227
        assert tree_counts=={'dead_all_prefixes':29226,'unknown_depth':1} and layers==266615 and maxdead==92 and maxlive==5
        # The preserved unknown is excluded analytically, never relabeled dead.
        from fractions import Fraction as F
        x=F(6,5);p=1+F(3,2)*x*x+3*x**3-F(3,2)*x**4-x**5+x**7;q=-(x*x+x**4)/2
        assert p>0 and q<0 and 13*q*q-p*p==F(446075579,6103515625)>0
    out={'status':'passed_full_S24_lattice_prefix_root_tree_and_positive_root_certificate','scope':'Degree7 real quadratic nonrational positions {2,4}, or reversal {3,5}; arbitrary finite signed integer cofactors. General conjecture open.','fields':53,'pair_count':774100,'row4_bases':33575,'row6_bases':20849,'prefix_counts':dict(prefix_counts),'root_counts':dict(root_counts),'tree_counts':dict(tree_counts),'tree_layers':layers,'deadmax':maxdead,'maxlive':maxlive,'unknown_excluded_by_positive_root':1,'workers':a.workers,'elapsed_seconds':time.monotonic()-started,'python':sys.version,'manifest_sha256':sha(r/'core-manifest.json'),'verifier_sha256':sha(Path(__file__))}
    out['direct_AB_replayed']=True
    partial=a.out.with_suffix(a.out.suffix+'.partial')
    with partial.open('w') as f:json.dump(out,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    partial.replace(a.out);print(json.dumps(out),flush=True)
if __name__=='__main__':main()
